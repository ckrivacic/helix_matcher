from helix.utils import numeric
from helix.utils.utils import max_subgraph
from helix.utils.utils import download_and_clean_pdb
import alphashape
from copy import deepcopy
import os
import numpy as np
import pandas as pd
import prody
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D
# from descartes import PolygonPatch


class ClashScore(object):
    def __init__(self, results_row, database_helices, query_helices,
            alpha=None):
        self.graph = results_row['graph']
        self.name = results_row['name'].split('_')[0]
        self.db_helices = database_helices
        self.query_helices = query_helices
        self.pdb_path = download_and_clean_pdb(self.name)
        print('PDB PATH for CLASH FILTER {}'.format(self.pdb_path))
        self.subgraphs = max_subgraph(self.graph)
        self.chain = self.db_helices.loc[self.subgraphs[0][0][0]]['chain']

        if not alpha:
            self.alpha = get_alphashape(self.pdb_path)
        else:
            self.alpha = alpha

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
            df_vectors, query_vectors = self.get_vectors(subgraph)
            transform = numeric.Transformation(df_vectors, query_vectors)
            prody_transform =\
                    prody.measure.transform.Transformation(transform.rotation,
                    transform.translation)
            prody_transform.apply(atoms)
            score = self.calculate(atoms)
            if score < best_score:
                best_score = score
                best_subgraph = subgraph

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

    def get_vectors(self, subgraph):
        '''Return the vectors for the helices in a subgraph'''
        df_vectors = []
        query_vectors = []
        for node in subgraph:
            df_idx = node[0]
            query_idx = node[1]
            df_row = self.db_helices.loc[df_idx]
            query_row = self.query_helices.loc[query_idx]

            df_vectors.extend(df_row['vector'])
            query_vectors.extend(query_row['vector'])    

        return df_vectors, query_vectors


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

    alpha_shape = alphashape.alphashape(coords, 0.18)

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
