from helix.utils import numeric
from helix.utils.utils import max_subgraph
from helix.utils.utils import download_and_clean_pdb
import alphashape
from copy import deepcopy
import math
import os
import numpy as np
import pandas as pd
import prody
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
# from descartes import PolygonPatch


def get_relative_path(workspace, path, depth=5):
    '''Take a path from the results dataframe (which is an absolute
    path) and get the path relative to the workspace root directory.
    This is a bit of a hack, and in the future, the paths should 
    probably be relative to begin with.'''
    # This should be everything after the root dir for a cluster
    # representative
    pathlist = path.split('/')[-depth:]
    pathlist.insert(0, workspace.root_dir)
    return os.path.join(*pathlist)


class ClashScore(object):
    def __init__(self, workspace, results_row, database_helices, query_helices,
            alpha=None, pdb=None, query_CAs=None):
        self.workspace = workspace
        self.graph = results_row['graph']
        self.name = results_row['name'].split('_')[0]
        self.db_helices = database_helices
        self.query_helices = query_helices
        if not pdb:
            self.pdb_path = download_and_clean_pdb(self.name)
        else:
            self.pdb_path = pdb
        print('PDB PATH for CLASH FILTER {}'.format(self.pdb_path))
        self.subgraphs = max_subgraph(self.graph)
        self.chain = self.db_helices.loc[self.subgraphs[0][0][0]]['chain']

        if not alpha:
            self.alpha = get_alphashape(self.pdb_path)
        else:
            self.alpha = alpha
        if not query_CAs:
            self.query_CAs = self.workspace.query_CAs

    def calc_interweave_score(self, atoms, df_rows, query_rows):
        '''
        Get indices of helices from df_rows to find which residues to
        look at.
        Load atoms for query helices form query_rows. (Maybe this should
        be in memory? Yeah. Load up query helices and pass them to this
        object.)
        Create a matrix of CA positions for scaffold atoms and query
        atoms. Or make a dataframe similar to patches, so we can sort?
        '''
        # tar_CAs = atoms.select('name CA').getCoords()

        def func(dist, w=20.0):
            return -w * math.cos(2 * dist * math.pi / 5.4) + w

        atoms = atoms.select("name CA")
        # score = 0
        scores = []
        for i in range(0, len(df_rows)):
            row = df_rows[i]
            target_helix_CAs = deepcopy(atoms)
            subselection = target_helix_CAs.select('resindex {}:{} and name '\
                    'CA'.format(row['start'] - 1, row['stop']))
            tar_CAs = subselection.getCoords()
            query_path = get_relative_path(self.workspace,
                    query_rows[i]['path'], depth=5)
            query_path = os.path.relpath(query_path,
                    start=self.workspace.root_dir)
            query_CAs = self.query_CAs[query_path]
            distance_matrix = pd.DataFrame()
            for i in range(0, len(query_CAs)):
                for j in range(0, len(tar_CAs)):
                    distance_matrix.at[i, j] = numeric.euclidean_distance(query_CAs[i], tar_CAs[j])
            helix_scores = []
            for j in range(0, len(tar_CAs)):
                nearest = distance_matrix[j].sort_values()
                nearest = nearest.iloc[0]
                if nearest < 5.4:
                    helix_scores.append(func(nearest))
            helix_score = sum(helix_scores) / len(helix_scores)
            scores.append(helix_score)

        return sum(scores)


    def apply(self):
        '''
        Find the score for each subgraph, return the score and subgraph
        of the best-scoring transformation
        '''
        best_score = 9999
        best_subgraph = None
        original_atoms = prody.parsePDB(self.pdb_path,
                chain=self.chain).select('backbone')
        for subgraph in self.subgraphs:
            atoms = deepcopy(original_atoms)
            # df_vectors, query_vectors = self.get_vectors(subgraph)
            df_rows, query_rows = self.get_helix_rows(subgraph)
            df_vectors = self.get_vectors(df_rows)
            query_vectors = self.get_vectors(query_rows)
            transform = numeric.Transformation(df_vectors, query_vectors)
            prody_transform =\
                    prody.measure.transform.Transformation(transform.rotation,
                    transform.translation)
            prody_transform.apply(atoms)
            score = self.calculate(atoms)
            if score < best_score:
                best_score = score
                best_subgraph = subgraph
            # self.interweave_score = self.calc_interweave_score(atoms, df_rows, query_rows)
            # assert(self.interweave_score is not None)
            # print(self.interweave_score)

        self.score = best_score
        self.subgraph = best_subgraph

    def calculate(self, atoms):
        '''Calculate a clash score from the transformed atoms object and
        alpha shape'''
        # For now just return the # of atoms inside the alpha shape
        clashes = 0
        for atom in atoms:
            coords = [atom.getCoords()]
            if self.alpha.contains(coords):
                clashes += 1

        return clashes

    def get_vectors(self, list_of_rows):
        vectors = []
        for row in list_of_rows:
            vectors.extend(row['vector'])
        return vectors


    def get_helix_rows(self, subgraph):
        '''Return the vectors for the helices in a subgraph'''
        # df_vectors = []
        # query_vectors = []
        df_rows = []
        query_rows = []
        for node in subgraph:
            df_idx = node[0]
            query_idx = node[1]
            df_rows.append(self.db_helices.loc[df_idx])
            query_rows.append(self.query_helices.loc[query_idx])
            # df_row = self.db_helices.loc[df_idx]
            # query_row = self.query_helices.loc[query_idx]

            # df_vectors.extend(df_row['vector'])
            # query_vectors.extend(query_row['vector'])    

        # return df_vectors, query_vectors
        return df_rows, query_rows


def get_alphashape(pdb, chain=None, plot=False):
    '''
    Returns an AlphaShape object of a pdb file, outlining its general
    shape for use in a clash filter.
    '''
    print('ALPHASHAPE PDB: {}'.format(pdb))
    if chain:
        atoms = prody.parsePDB(pdb, chain=chain)
    else:
        atoms = prody.parsePDB(pdb)
        atoms = atoms.select('not chain A')
    # For some reason there is a level which is not populated. This may
    # be for proteins w/ multiple chains.
    coordsets = atoms.getCoordsets()
    coords = []
    for coordset in coordsets:
        coords.extend(coordset)

    # coords = [(0., 0.), (0., 1.), (1., 1.), (1., 0.), (0.5, 0.5)]

    alpha_shape = alphashape.alphashape(coords, 0.1)

    if plot:
        helix = prody.parsePDB(pdb, chain='A')
        helixcoords = helix.getCoordsets()[0]
        fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        ax = Axes3D(fig)
        ax.scatter(*zip(*coords))
        ax.scatter(*zip(*helixcoords))
        # # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
        ax.plot_trisurf(*zip(*alpha_shape.vertices),
                triangles=alpha_shape.faces, alpha=0.3)
        plt.show()
    return alpha_shape


def main():
    pdb = os.path.abspath('test_files/boundary/cluster_representatives/clst10_rep.pdb.gz')
    alpha = get_alphashape(pdb, chain='B')

def test():
    helixdir = os.path.expanduser('~/intelligent_design/helix_matcher/')
    query_helices = os.path.join(helixdir, 'test_files/boundary/cluster_representatives/query_helices.pkl')
    query_helices = pd.read_pickle(query_helices)
    df_helices = os.path.join(helixdir, 'dataframes_clash/final.pkl')
    df_helices = pd.read_pickle(df_helices)
    results = os.path.join(helixdir, 'clashout_30/query_results_000.pkl')
    results = pd.read_pickle(results)
    results.sort_values(by='matches', inplace=True, ascending=False)

    pdb = os.path.join(helixdir, 'test_files/boundary/cluster_representatives/clst10_rep.pdb.gz')
    alpha_shape = get_alphashape(pdb, chain=None, plot=True)

    for idx, row in results.iterrows():
        clash_score = ClashScore(row, df_helices, query_helices,
                alpha=alpha_shape)
        clash_score.apply()
        print(row['name'])
        print(clash_score.score)


if __name__=='__main__':
    test()
# main()
