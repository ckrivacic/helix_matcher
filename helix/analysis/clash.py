from helix.utils import numeric
from helix.utils.utils import max_subgraph
from helix.utils.utils import download_and_clean_pdb
from helix.utils.geometry import vector_intersects_plane
import alphashape
from copy import deepcopy
import math
import os
import numpy as np
import pandas as pd
import prody
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
from mpl_toolkits.mplot3d import axes3d, Axes3D
# from descartes import PolygonPatch


def get_relative_path(workspace, path, depth=5):
    '''Take a path from the results dataframe (which is an absolute
    path) and get the path relative to the workspace root directory.
    This is a bit of a hack, and in the future, the paths should 
    probably be relative to begin with.'''
    # This should be everything after the workspace root dir for a cluster
    # representative
    pathlist = path.split('/')[-depth:]
    pathlist.insert(0, workspace.root_dir)
    return os.path.join(*pathlist)


class Score(object):
    def __init__(self, workspace, results_row, database_helices, query_helices,
            alpha=None, query_CAs=None, target_path=None):
        self.workspace = workspace
        # Graph object indicating which indices are matched
        self.graph = results_row['graph']
        self.name = results_row['name']#.split('_')[0]
        # Dataframe of matched helices
        self.db_helices = database_helices
        # All dense subgraphs that have the maximum number of members
        self.subgraphs = max_subgraph(self.graph)
        # Dataframe of all docked helices
        self.query_helices = query_helices
        if 'path' not in database_helices.columns:
            print(database_helices.columns)
            self.pdb_path = download_and_clean_pdb(self.name.split('_')[0])
        else:
            self.pdb_path = os.path.join(
                    self.workspace.root_dir,
                    self.db_helices.loc[self.subgraphs[0][0][0]]['path']
                    )
        print('PDB PATH for CLASH FILTER {}'.format(self.pdb_path))
        # Chain of database helices
        self.chain = self.db_helices.loc[self.subgraphs[0][0][0]]['chain']

        # Make an alphashape if one is not provided
        if not alpha:
            if not target_path:
                print('No alphashape or target provided. Do not attempt to '\
                        'calculate clash score.')
            else:
                self.alpha = get_alphashape(self.target_path)
        else:
            self.alpha = alpha
        if not query_CAs:
            # PRODY CAs of all query helices
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

        def func(dist, w=15):
            return -w * math.cos(2 * dist * math.pi / 5.4) + w

        atoms = atoms.select("name CA")
        # score = 0
        scores = []
        for i in range(0, len(df_rows)):
            row = df_rows[i]
            target_helix_CAs = atoms.copy()
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
                if nearest < 5.4:# and nearest > 1.3:
                    # helix_scores.append(nearest)
                    helix_scores.append(func(nearest))
                # else:
                    # helix_scores.append(0)
            if len(helix_scores) > 0:
                helix_score = sum(helix_scores) / len(helix_scores)
                scores.append(helix_score)
            else:
                scores.append(-1)

        return sum(scores)

    def parallel_rmsd(self, atoms, df_rows, query_rows, transformer):
        '''Finds the RMSD of only the overlapping portions of helices'''
        rmsds = []
        lengths = []
        for i in range(0, len(df_rows)):
            row = df_rows[i]
            vector = row['vector']
            vector = numeric.apply_transformation(transformer, vector)
            target_helix_CAs = atoms.copy()
            subselection = target_helix_CAs.select('resindex {}:{} and name' \
                                                   ' CA'.format(row['start'] - 1, row['stop']))
            tar_CAs = subselection.getCoords()

            query_path = get_relative_path(self.workspace,
                                           query_rows[i]['path'], depth=5)
            query_path = os.path.relpath(query_path,
                                         start=self.workspace.root_dir)
            query_CAs = self.query_CAs[query_path]

            start_ca_query = query_CAs[0]
            overlap_query = []
            for atom in range(1, len(query_CAs)):
                query_vector = np.array([start_ca_query, query_CAs[atom]])
                # NOTE: Need to add logic to handle "none" returns, or change the function to return >1 if parallel.
                if vector_intersects_plane(query_vector, vector[0], vector) < 1 and \
                        vector_intersects_plane(query_vector, vector[1], vector) > 1:
                    overlap_query.append(query_CAs[atom])
            overlap_scaffold = []
            start_ca_scaffold = tar_CAs[0]
            for atom in range(1, len(tar_CAs)):
                scaffold_vector = np.array([start_ca_scaffold, tar_CAs[atom]])
                # NOTE: Need to add logic to handle "none" returns, or change the function to return >1 if parallel.
                if vector_intersects_plane(scaffold_vector, query_vector[0], query_vector) < 1 and \
                        vector_intersects_plane(scaffold_vector, query_vector[1], query_vector) > 1:
                    overlap_scaffold.append(tar_CAs[atom])

            rmsd = find_best_rmsd(overlap_query, overlap_scaffold)
            if len(overlap_scaffold) > len(overlap_query):
                length = len(overlap_query)
            else:
                length = len(overlap_scaffold)
            rmsds.append(rmsd)
            lengths.append(length)

        weighted = 0
        for rmsd, length in zip(rmsds, lengths):
            weighted += rmsd * length

        return weighted / sum(lengths)

    def calc_rmsd(self, atoms, df_rows, query_rows):
        '''Calculate the best possible RMSD of each matched helix to the docked helices'''
        # Requires database atoms, row, and query atoms objects
        # "atoms" will be database atoms
        lengths = []
        rmsds = []
        for i in range(0, len(df_rows)):
            row = df_rows[i]
            target_helix_CAs = atoms.copy()
            subselection = target_helix_CAs.select('resindex {}:{} and name'\
                                                   ' CA'.format(row['start'] - 1, row['stop']))
            tar_CAs = subselection.getCoords()
            query_path = get_relative_path(self.workspace,
                                           query_rows[i]['path'], depth=5)
            query_path = os.path.relpath(query_path,
                                         start=self.workspace.root_dir)
            query_CAs = self.query_CAs[query_path]

            rmsd = find_best_rmsd(tar_CAs, query_CAs)
            if len(tar_CAs) > len(query_CAs):
                length = len(query_CAs)
            else:
                length = len(tar_CAs)
            rmsds.append(rmsd)
            lengths.append(length)

        weighted = 0
        for rmsd, length in zip(rmsds, lengths):
            weighted += rmsd * length

        return weighted / sum(lengths)

    def apply(self):
        '''
        Find the score for each subgraph, return the score and subgraph
        of the best-scoring transformation
        '''
        best_score = 9999
        best_subgraph = None
        print(f'OPENING PDBPATH {self.pdb_path}')
        original_atoms = prody.parsePDB(self.pdb_path,
                chain=self.chain).select('backbone')
        for subgraph in self.subgraphs:
            atoms = original_atoms.copy()
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
            interweave_score = self.calc_interweave_score(atoms, df_rows, query_rows)
            # assert(self.interweave_score is not None)
            # print(self.interweave_score)
            if score + interweave_score < best_score:
                best_interweave_score = interweave_score
                best_clash_score = score
                best_score = score + interweave_score
                best_subgraph = subgraph

        self.interweave_score = best_interweave_score
        self.score = best_clash_score
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


def find_best_rmsd(CAs_1, CAs_2):
    '''Calculate the best possible RMSD between two sets of atoms (intended to be CAs only)
    without moving either of them'''
    if len(CAs_2) > len(CAs_1):
        longer = CAs_2
        shorter = CAs_1
        flipped = True
    else:
        longer = CAs_1
        shorter = CAs_2
        flipped = False

    long_combos = []
    for i in range(1, len(longer) - len(shorter) + 2):
        # Say longer is 5 atoms, shorter is 4 atoms; 5 - 4 = 1, but there are 2 possibilities for overlap
        resi_stop = i + len(shorter) - 1
        # long_combos.append(longer.select(f'resindex {i} to {resi_stop}'))
        long_combos.append(longer[i-1:resi_stop])
    short_combos = [shorter]

    best_rmsd = 99999
    # short_seq = None
    # long_seq = None
    # short_resis = (-1, -1)
    # long_resis = (-1, -1)
    for long_sele in long_combos:
        for short_sele in short_combos:
            if len(long_sele) == len(short_sele):
                # rmsd = prody.calcRMSD(long_sele, short_sele)
                rmsd = numeric.RMSD(long_sele, short_sele)
            else:
                print("\033[91mError: atom lengths did not match.\033[0m")
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                # short_seq = short_sele.getSequence()
                # long_seq = long_sele.getSequence()
    if best_rmsd == 99999:
        print('Error with best RMSD')
        print('Short combos:', flush=True)
        print(long_combos, flush=True)
        print('Long combos:', flush=True)
        print(short_combos, flush=True)
        print('CAs_1:', flush=True)
        print(CAs_1, flush=True)
        print('CAs_2:', flush=True)
        print(CAs_2, flush=True)

    return best_rmsd


def get_alphashape(pdb, chain=None, plot=False):
    '''
    Returns an AlphaShape object of a pdb file, outlining its general
    shape for use in a clash filter.
    '''
    warnings.filterwarnings("ignore", message="Mesh is non-watertight for contained point query!")

    print('ALPHASHAPE PDB: {}'.format(pdb))
    if chain:
        atoms = prody.parsePDB(pdb, chain=chain)
    else:
        atoms = prody.parsePDB(pdb)
        atoms = atoms.select('chain A')
    # For some reason there is a level which is not populated. This may
    # be for proteins w/ multiple chains.
    coordsets = atoms.getCoordsets()
    coords = []
    for coordset in coordsets:
        coords.extend(coordset)

    # coords = [(0., 0.), (0., 1.), (1., 1.), (1., 0.), (0.5, 0.5)]

    alpha = 0.2
    alpha_shape = alphashape.alphashape(coords, alpha)

    if plot:
        mpl.use('tkagg')
        helix = prody.parsePDB(pdb, chain='A')
        helixcoords = helix.getCoordsets()[0]
        fig = plt.figure()
        # ax = fig.add_subplot(projection='3d')
        ax = Axes3D(fig)
        ax.scatter(*zip(*coords))
        ax.scatter(*zip(*helixcoords))
        # # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2))
        ax.plot_trisurf(*zip(*alpha_shape.vertices),
                triangles=alpha_shape.faces, alpha=0.4)
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
        clash_score = Score(row, df_helices, query_helices,
                alpha=alpha_shape)
        clash_score.apply()
        print(row['name'])
        print(clash_score.score)


def test_rmsd():
    import helix.workspace as ws
    from helix.utils.utils import safe_load
    workspace = ws.workspace_from_dir(os.path.expanduser('~/intelligent_design/2022_01_05_denovo/rifdock_outputs/5U8R/'))
    query_helices = os.path.join(workspace.focus_dir, 'filtered', 'query_helices.pkl')
    query_helices = safe_load(query_helices)
    df_helices = os.path.join(workspace.project_params_dir, 'database', 'helixdf_custom_08-13-2021.pkl')
    df_helices = safe_load(df_helices)


if __name__=='__main__':
    test()
# main()
