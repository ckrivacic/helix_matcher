'''
Create bins or match a query protein.
Usage:
    matcher.py bin <helix_dataframe> [options]
    matcher.py match <pdb> [options]

options:
    --local, -l  Run locally
    --verbose, -v  Verbose output
'''
import docopt
import pandas as pd
import numpy as np
import numeric
from itertools import product
import os, psutil, sys
import pickle
import subprocess
from scan_helices import final_vector
from pymongo import MongoClient
from pyrosetta import init, pose_from_file
# import networkx as nx
import graph_tool.all as gt
import collections


def plot_vectors(vectors, color='darkgray'):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for vector in vectors:
        x = [point[0] for point in vector]
        y = [point[1] for point in vector]
        z = [point[2] for point in vector]
        ax.plot(x, y, z, color=color, linewidth=4)
    plt.show()


def bin_array(array, bins):
    '''
    Digitize a numpy array.
    TO DO: Circularize the binning of angles somehow.
    '''
    inds = np.digitize(array, bins)
    binned = tuple([bins[inds[n]-1] for n in range(array.size)])
    return binned


def relative_position(row1, row2, vectortype='normalized_vector'):
    '''
    Gives the internal relative orientation of two lines, given their
    row from the pandas dataframe created in scan_helices.
    The relative orientation of two lines should be able to be described
    with just 4 parameters, since they are 2D objects in 3D space. If we
    have lines consisting of points [a,b] and [c,d], those parameters are:
    - The distance between their centroids
    - Angle abc
    - Angle bcd
    - Dihedral abcd
    '''

    norm_v1 = row1[vectortype]
    norm_v2 = row2[vectortype]
    centroid_dist = numeric.euclidean_distance(row1['centroid'],
            row2['centroid'])
    abc = numeric.angle(norm_v1[0], norm_v1[1], norm_v2[0])
    bcd = numeric.angle(norm_v1[1], norm_v2[0], norm_v2[1])
    dihedral = numeric.dihedral(norm_v1[0], norm_v1[1], norm_v2[0],
            norm_v2[1])
    # plot_vectors([norm_v1, norm_v2], color='black')

    return centroid_dist, abc, bcd, dihedral


class Match(object):
    '''
    Class to construct a potential match.
    '''
    def __init__(self, name, query_db, main_db):
        self.name = name
        self.query = query_db
        self.db = main_db.xs(name, level='name')
        self.graph = gt.Graph(directed=False)
        # Track helix pairs so we don't add them to the graph more than
        # once

    def max_subgraph(self):
        '''
        Finds dense subgraphs, which represent compatible sets of helix
        pairs between the query helices and the database PDB. The
        longest such subgraph represents the best overlay of the PDB
        with the set of query helices.
        '''
        max_subgraph_len = 0
        for f in gt.max_cliques(self.graph):
            if len(f) > max_subgraph_len:
                max_subgraph_len = len(f)

        print('Max number of matches:')
        print(max_subgraph_len)
        return max_subgraph_len


    def plot_graph(self):
        # import matplotlib.pyplot as plt
        # import graph_tool.draw as draw
        # plt.subplot(111)
        gt.remove_parallel_edges(self.graph)
        pos = gt.fruchterman_reingold_layout(self.graph, n_iter=1000)
        gt.graph_draw(self.graph, pos=pos)
        # plt.show()

    def find_edges(self):
        '''
        Populate the graph with nodes and edges.
        Each node consists of a pair of indices, one from the main
        database and one from the query database. This pairing
        represents the case where the helix in the first index is
        overlaid on the helix of the second index. Edges represent
        compatibility between adjacent nodes.
        '''
        print('Finding edges')
        edges = []
        self.nodes = {}
        i = 0
        for doc in self.db.iterrows():
            compatible_bins = self.query.xs(doc.index())
            # compatible_bins = self.query.find({'bin': doc['bin']})
            for result in compatible_bins.iterrows():
                idx_pair1 = (doc[1]['idx1'], result[1]['idx1'])
                idx_pair2 = (doc[1]['idx2'], result[1]['idx2'])
                # Track which nodes have been sampled
                if idx_pair1 not in self.nodes:
                    self.nodes[idx_pair1] = i
                    i += 1
                    # self.nodes.append(idx_pair1)
                    # self.graph.add_node(idx_pair1)
                if idx_pair2 not in self.nodes:
                    self.nodes[idx_pair2] = i
                    i += 1
                    # self.nodes.append(idx_pair2)
                    # self.graph.add_node(idx_pair2)
                # print('Edge found:')
                # print(idx_pair1)
                # print(idx_pair2)
                edges.append((self.nodes[idx_pair1],
                    self.nodes[idx_pair2]))
                # i += 2
        # nodes = set(self.nodes)
        # self.graph.add_edge(idx_pair1, idx_pair2)
        # print(nodes)
        print('All edges:')
        print(edges)
        self.graph.add_edge_list(edges)


class HelixBin(object):
    def __init__(self, helix_db, exposed_cutoff=0.3, length_cutoff=10.8,
            query_df=None, query_name=None, angstroms=2.5, degrees=15,
            verbose=False, start=None, stop=None):
        self.verbose = verbose
        self.df = helix_db
        self.df['idx'] = self.df.index

        # Binning parameters
        self.degrees = degrees
        self.angstroms = angstroms
        self.setup_bins()
        binned_name = 'bins_{}A_{}D'.format(self.angstroms,
                self.degrees)
        self.start = start
        self.stop = stop

        # Trimming dataframe
        if length_cutoff:
            self.df = self.df[self.df['length'] > length_cutoff]
        if exposed_cutoff:
            self.df = self.df[self.df['percent_exposed'] >
                    exposed_cutoff]
        if 'normalized_vector' not in self.df.columns:
            self.df['normalized_vector'] = self.df.apply(lambda x:
                    final_vector(x['direction'], 1, x['centroid']), axis=1)

    def setup_bins(self):
        nrbins = int(360//self.degrees) + 1
        self.rbins = np.linspace(-180, 180, nrbins)
        tstart = -10000
        tstop = 10000
        ntbins = int((tstop - tstart) // self.angstroms) + 1
        self.tbins = np.linspace(tstart, tstop, ntbins)

    def bin_db(self, check_dups=False):
        '''
        Bin dataframes.
        '''

        from scipy.spatial.transform import Rotation as R
        import subprocess
        import time

        # db = self.client[dbname]
        # bins = db['bins_{}A_{}D'.format(
            # self.angstroms, self.degrees
            # )]
        bins = pd.DataFrame(columns=['bin', 'name', 'idx1', 'idx2'])
        # Pandas indices are hash lookups and we can have multiple of
        # them, but they cannot be added piecewise. Therefore we will
        # create partial tables, then create the indices and save the
        # dataframes. Results will be saved in chunks.
        # bins.set_index(['bin', 'name'], inplace=True)
        total_proteins = len(set(self.df['name']))
        interval = 500

        # import shelve

        # binned = shelve.open('binned_0p3/hashtable', 'c', writeback=True)
        # i tracks # of names analyzed
        i = 0
        # saveno tracks how many dataframes have been saved.
        self.saveno = 1
        unsaved_docs = []
        start_time = time.time()

        def update(bins, start_time, unsaved_docs, interval, i,
                final=False):
            print('{} of {} PDBs processed so far.'.format(
                i, total_proteins))
            mem_used = psutil.Process(os.getpid()).memory_info().rss
            if self.verbose:
                print('Currently using {} GB of memory'.format(
                    mem_used * 10**-9
                    ))
            df_mem = bins.memory_usage(index=True, deep=True).sum()
            if self.verbose:
                print('Dataframe is using {} GB of memory'.format(
                    df_mem * 10**-9
                    ))
            elapsed = time.time() - start_time
            rate = interval / elapsed
            remaining = (total_proteins - i) / rate / 3600
            print('Analysis of 500 pdbs took {} seconds. Est. {} h remaining'.format(
                elapsed, remaining
                ))

            if len(unsaved_docs) > 0:
                if self.verbose:
                    print('Adding to dataframe...')
                bins = bins.append(unsaved_docs, ignore_index=True)
                if self.verbose:
                    print(bins)
            else:
                if self.verbose:
                    print('Nothing to update for this batch.')


            # Save when memory footprint of dataframe gets larger than 4
            # GB. This way each sub-dataframe can be read into memory.
            if df_mem * 10**-9 > 4 or final:
                bins.set_index(['bin', 'name'], inplace=True)
                outfolder = 'database/bins_{}A_{}D/'.format(self.angstroms, self.degrees)
                outfile = 'bins_{}A_{}D_{:04d}.pkl'.format(self.angstroms,
                        self.degrees, self.saveno)
                out = os.path.join(outfolder, outfile)
                print('Saving current dataframe to {}'.format(out))
                if not os.path.exists(outfolder):
                    os.makedirs(outfolder, exist_ok=True)
                bins.to_pickle(out)
                self.saveno += 1
                if self.verbose:
                    print('Saved.')

                # If saved to disk, return an empty dataframe.
                return pd.DataFrame()

            else:
                # Return input dataframe if we have not saved it to disk.
                return bins

        groups = self.df.groupby(['name'])
        names = sorted(list(groups.groups.keys()))
        if self.start:
            names = names[self.start:]
        if self.stop:
            names = names[:self.stop]
        for name in names:
        # for name, group in df.groupby(['name']):
            group = groups.groups[name]
            i += 1

            for combination in product(self.df.loc[group].T.to_dict().values(),
                    repeat=2):
                if combination[0]['idx'] != combination[1]['idx']:
                    # vector1 = combination[0]['vector']
                    # vector2 = combination[1]['vector']

                    # plot_vectors([vector1, vector2], color='purple')

                    idx1 = combination[0]['idx']
                    idx2 = combination[1]['idx']
                    # if self.verbose:
                        # print('------------------------------------')
                        # print(combination[0])
                        # print(combination[1])

                    dist, angle1, angle2, dihedral =\
                            relative_position(combination[0], combination[1])
                    dist = np.array([dist])
                    angles = np.array([angle1, angle2, dihedral])

                    rbin = bin_array(angles, self.rbins)
                    tbin = bin_array(dist, self.tbins)
                    rbin2 = bin_array(angles, self.rbins + (self.degrees/2))
                    tbin2 = bin_array(dist, self.tbins +
                            (self.angstroms/2))

                    x = [tbin[0], tbin2[0]]
                    abc = [rbin[0], rbin2[0]]
                    bcd = [rbin[1], rbin2[1]]
                    dih = [rbin[2], rbin2[2]]

                    for bin_12 in product(x, abc, bcd,
                        dih):
                        bin_12 = ' '.join(map(str, bin_12))
                        doc = {
                                'bin':bin_12, 
                                'name': name,
                                'idx1':idx1,
                                'idx2':idx2
                        }
                        # if check_dups:
                            # if len(list(bins.find(doc))) == 0:
                                # unsaved_docs.append(doc)
                        # else:
                        unsaved_docs.append(doc)

            if i%interval == 0:
                bins = update(bins, start_time, unsaved_docs, interval, i)
                start_time = time.time()
                unsaved_docs = []

        bins = update(bins, start_time, unsaved_docs, interval, i, final=True)

        return bins

class HelixLookup(object):
    '''
    Class to handle binning and matching of helix databases. This maybe
    should be two classes, one for binning and one for matching, but
    this is it for now.
    '''

    def __init__(self, lookup, query):
        self.lookup = lookup
        self.query = query

    def score_match(self, list_of_index_pairs):
        """
        Idea (idk where else to put this):
            To get 3rd, 4th, etc. helices, do a reverse lookup. That is,
            for each bin in the FOUND PDB, look for matches in the QUERY
            pdb.
        """
        for tup in list_of_index_pairs:
            query_row = self.query_df.loc[tup[0]]
            db_row = self.df.loc[tupe[1]]


    def match(self):
        names = []

        # Pandas rewrite
        for _bin, group in self.query.groupby(level='bin'):
            for result in self.lookup.xs(_bin, level='bin').iterrows():
                names.append(
                        result['name']
                        )

        print('Forward search done.')

        min_matches = 4
        names = [item for item, count in
                collections.Counter(names).items() if
                count >= min_matches]

        print(names)
        print(len(names))

        results = []
        # TEMP

        # sys.exit()
        i = 0
        for name in names:
            i += 1
            result = {}
            result['name'] = name
            print('-------------------------------------------------')
            print('Name: {}'.format(name))
            match = Match(name, self.query, self.lookup)
            match.find_edges()
            result['matches'] = match.max_subgraph()
            result['graph'] = match.graph
            results.append(result)
            if i % 100 == 0:
                df = pd.DataFrame(results[i-99:i+1])
                df.to_pickle('match_results_{}.pkl'.format(i))
            # match.plot_graph()
            # print('searching {}'.format(name))
            # for _bin in self.binned.find({'name': name[0]}):
                # if _bin['idx1'] == name[1]:
                    # print('-------')
                    # print(_bin)
                    # for doc in self.query_bins.find({'bin':_bin['bin']}):
                        # print('MATCH:')
                        # results[name].append((doc['idx1'], doc['idx2']))
                        # print(doc)

        df = pd.DataFrame(results)
        df.to_pickle('match_results.pkl')
        # for key in results:
            # print('------------------RESULTS FOR {}----------------'.format(
                            # key
                        # ))
            # for pair in set(results[key]):
                # print(pair)
        # for key in results:
            # print('PDB {} had {} matching transformations'.format(
                # key, len(set(results[key]))
                # ))


def test():
    import scan_helices

    test_path = 'test_files/6r9d.cif'
    init()
    pose = pose_from_file(test_path).split_by_chain(2)
    print(pose.size())
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices()
    helices = pd.DataFrame(helices)
    print(helices)
    helices = helices[helices['percent_exposed'] > 0.3]
    print(helices)
    print(helices.shape)
    print(helices['name'])

    # lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            # query_df=helices, query_name='6r9d')
    lookup = HelixLookup(pd.DataFrame(),
            query_df=helices, query_name='6r9d', angstroms=5,
            # degrees=15, reset_querydb=True, dbname='nr')
            degrees=30, reset_querydb=True, dbname='test_bins')
    lookup.match()

def test_rifdock():
    import scan_helices

    test_path = 'test_files/test_rifgen/cluster_representatives/matchme.pdb'
    init()
    pose = pose_from_file(test_path)
    print(pose.size())
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices(split_chains=False,
            name='rifdock_test')
    helices = pd.DataFrame(helices)
    helices.to_pickle('rifdock_helices.pkl')
    sys.exit()
    print(helices)
    # helices = helices[helices['percent_exposed'] > 0.3]
    print(helices)
    print(helices.shape)
    print(helices['name'])

    # lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            # query_df=helices, query_name='6r9d')
    lookup = HelixLookup(pd.DataFrame(),
            query_df=helices, query_name='6r9d', angstroms=2.5,
            degrees=15, reset_querydb=True, dbname='nr')
            # degrees=30, reset_querydb=True, dbname='test_bins')
    lookup.match()

def make_hash_table():
    print('Loading database and setting up lookup object...')
    # length cutoff of 2 turns or 10.8 angstroms
    lookup = HelixLookup(pd.read_pickle('nr_dataframes/final.pkl'),
            exposed_cutoff=0.3, length_cutoff=10.8, angstroms=2.5,
            degrees=15, dbname='nr')
    print('Done.')
    # binned = lookup.bin_db(lookup.df)
    lookup.update_bin_db()
    # out = "binned_0p3/last.pkl"
    # with open(out, 'wb') as f:
        # pickle.dump(binned, f)

def make_test_hash_table():
    client = MongoClient()
    deg=15
    angstroms=2.5
    # client['test_bins']['bins_{}A_{}D'.format(angstroms, deg)].drop()
    lookup=HelixLookup(pd.read_pickle('out.pkl'), exposed_cutoff=0.3,
            length_cutoff=10.8, angstroms=angstroms, degrees=deg,
            dbname='test_bins')
    lookup.update_bin_db()


def main():
    args = docopt.docopt(__doc__)
    if args['bin']:
        lookup = HelixBin(pd.read_pickle(args['<helix_dataframe>']),
                exposed_cutoff=0.3, length_cutoff=10.8, angstroms=2.5,
                degrees=15, verbose=args['--verbose'])
        lookup.bin_db()




if __name__=='__main__':
    # test()
    # test_rifdock()
    # make_hash_table()
    # make_test_hash_table()
    main()
