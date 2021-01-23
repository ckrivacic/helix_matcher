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
import networkx as nx
'''
Here's the plan.

1. Hash the table? Idk. Look into this.
2. Take a pair of helices and align it to a pair of helices in the
database.
3. Look for 3rd, 4th matching helices.

Okay, how to see if a helix matches?
Both tips within a certain distance? EZPZ.


Hashing:
    1. For PDB database, take all pairs of helices within a protein and
    bin their relative positions. (Need to bin twice for each
    orientation and distance.)
    2. Hash binned relative positions and store the indices from the
    original PDB database. Maybe limit hashing to helices with some %
    surface exposure.
    3. Hash query helix pair, lookup hash table, then lookup original
    helices.

A dictionary is already a hash table. So just use the binned
superposition matrices as the keys.

'''


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
        self.db = main_db.find({'name':name})
        self.graph = nx.Graph()
        # Track helix pairs so we don't add them to the graph more than
        # once
        self.seen_nodes = set()

    def max_subgraph(self):
        '''
        Finds dense subgraphs, which represent compatible sets of helix
        pairs between the query helices and the database PDB. The
        longest such subgraph represents the best overlay of the PDB
        with the set of query helices.
        '''
        for f in nx.find_cliques(self.graph):
            print(f)

    def plot_graph(self):
        import matplotlib.pyplot as plt
        plt.subplot(111)
        nx.draw(self.graph, with_labels=True, font_weight='bold')
        plt.show()

    def find_edges(self):
        '''
        Populate the graph with nodes and edges.
        Each node consists of a pair of indices, one from the main
        database and one from the query database. This pairing
        represents the case where the helix in the first index is
        overlaid on the helix of the second index. Edges represent
        compatibility between adjacent nodes.
        '''
        for doc in self.db:
            compatible_bins = self.query.find({'bin': doc['bin']})
            for result in compatible_bins:
                idx_pair1 = (doc['idx1'], result['idx1'])
                idx_pair2 = (doc['idx2'], result['idx2'])
                # Track which nodes have been sampled
                if idx_pair1 not in self.seen_nodes:
                    self.seen_nodes.add(idx_pair1)
                    self.graph.add_node(idx_pair1)
                if idx_pair2 not in self.seen_nodes:
                    self.seen_nodes.add(idx_pair2)
                    self.graph.add_node(idx_pair2)
                self.graph.add_edge(idx_pair1, idx_pair2)



class HelixLookup(object):
    '''
    Class to handle binning and matching of helix databases. This maybe
    should be two classes, one for binning and one for matching, but
    this is it for now.
    '''

    def __init__(self, df, exposed_cutoff=0.5, length_cutoff=10.8,
            query_df=None, query_name=None, angstroms=2, degrees=20,
            reset_querydb=False, dbname='helix_bins', verbose=False):
        # Setup pymongo
        self.client = MongoClient()
        try:
            print(self.client.server_info())
        except:
            print('Problem with MongoDB server connection; attempting to open...')
            subprocess.Popen(['mongod', '--dbpath',
                '/Users/codykrivacic/data/db'])
            print(self.client.server_info())

        self.verbose = verbose
        # Dataframe consisting of vectors (is this needed when not
        # creating database?)
        self.df = df
        self.df['idx'] = self.df.index
        if 'normalized_vector' not in self.df.columns:
            self.df['normalized_vector'] = self.df.apply(lambda x:
                    final_vector(x['direction'], 1, x['centroid']), axis=1)
        # TEMPORARY try/except just to speed things up (not feed a
        # dataframe)
        try:
            if exposed_cutoff is not None:
                self.df = self.df[self.df['percent_exposed'] > exposed_cutoff]
            if length_cutoff is not None:
                self.df = self.df[self.df['length'] > length_cutoff]
        except:
            pass
        self.degrees = degrees
        self.angstroms = angstroms
        self.setup_bins()
        binned_name = 'bins_{}A_{}D'.format(self.angstroms, self.degrees)

        # Questioning whether to define this here or not.
        self.dbname = dbname
        self.binned = self.client[self.dbname][binned_name]

        # Setting up the query df, if provided.
        if query_df is not None:
            self.query_df = query_df
            if query_name is None:
                query_name = self.query_df.iloc[0]['name']
            self.query_df['idx'] = self.query_df.index
            if 'normalized_vector' not in self.query_df.columns:
                self.query_df['normalized_vector'] =\
                        self.query_df.apply(lambda x:
                            final_vector(x['direction'], 1, x['centroid']),
                            axis=1)
            if reset_querydb:
                print('Deleting old query db')
                client = MongoClient()
                client[query_name][binned_name].drop()
            self.query_bins = self.bin_db(self.query_df, query_name,
                    check_dups=True)

    def setup_bins(self):
        nrbins = int(360//self.degrees) + 1
        self.rbins = np.linspace(-180, 180, nrbins)
        tstart = -10000
        tstop = 10000
        ntbins = int((tstop - tstart) // self.angstroms) + 1
        self.tbins = np.linspace(tstart, tstop, ntbins)

    def update_bin_db(self):
        self.bin_db(self.df, self.dbname, check_name=True)

    def bin_db(self, df, dbname, check_name=True, check_dups=False):
        '''
        Save parameter tells pymongo where to save (must have MongoDB
        installed and running).
        '''

        from scipy.spatial.transform import Rotation as R
        import subprocess
        import time

        db = self.client[dbname]
        bins = db['bins_{}A_{}D'.format(
            self.angstroms, self.degrees
            )]
        total_proteins = len(set(df['name']))
        interval = 500

        # import shelve

        # binned = shelve.open('binned_0p3/hashtable', 'c', writeback=True)
        i = 0
        # j = 0
        unsaved_docs = []
        start_time = time.time()

        def update(bins, start_time, unsaved_docs, interval, i):
            mem_used = psutil.Process(os.getpid()).memory_info().rss
            print('{} of {} PDBs processed so far.'.format(
                i, total_proteins))
            print('Currently using {} G of memory'.format(
                mem_used * 10**-9
                ))
            print('Saving to db...')
            if len(unsaved_docs) > 0:
                bins.insert_many(unsaved_docs, ordered=False)
                bins.create_index([('bin', 'hashed')])
                bins.create_index([('name', 'hashed')])
            else:
                print('Nothing to update for this batch.')
            elapsed = time.time() - start_time
            rate = interval / elapsed
            remaining = (total_proteins - i) / rate / 3600
            print('Done. 500 pdbs took {} seconds. Est. {} h remaining'.format(
                elapsed, remaining
                ))

        for name, group in df.groupby(['name']):
            if check_name:
                if len(list(bins.find({'name':name}))) > 0:
                    if i%interval == 0:
                        i += 1
                        update(bins, start_time, unsaved_docs, interval, i)
                        start_time = time.time()
                        unsaved_docs = []
                    continue

            i += 1
            for combination in product(group.T.to_dict().values(),
                    repeat=2):
                if combination[0] != combination[1]:
                    # vector1 = combination[0]['vector']
                    # vector2 = combination[1]['vector']

                    # plot_vectors([vector1, vector2], color='purple')

                    idx1 = combination[0]['idx']
                    idx2 = combination[1]['idx']
                    # transform = numeric.Transformation(vector1,
                            # vector2)
                    # rot = R.from_matrix(transform.rotation)
                    if self.verbose:
                        print('------------------------------------')
                        print(combination[0])
                        print(combination[1])
                    # rot = rot.as_euler('ZYZ', degrees=True)
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
                        if check_dups:
                            if len(list(bins.find(doc))) == 0:
                                unsaved_docs.append(doc)
                        else:
                            unsaved_docs.append(doc)

                        # bins.insert_one(doc)
                        # if bin_12 not in binned:
                            # binned[bin_12] = {}
                        # if name not in binned[bin_12]:
                            # binned[bin_12][name] = []
                        # binned[bin_12][name].append((idx1, idx2))
            if i%interval == 0:
                update(bins, start_time, unsaved_docs, interval, i)
                start_time = time.time()
                unsaved_docs = []
                # For large databases, dump dictionaries when memory
                # usage goes over 6 GB. Only use for forming initial lookup
                # database.
                # print('Syncing db...')
                # binned.sync()
                # print('Done.')
                # If more than 6G memory used, dump just in case.
                # if mem_used > 12 * 10**9:
                    # out = "binned_0p3/{}.pkl".format(j)
                    # j += 1
                    # with open(out, 'wb') as f:
                        # pickle.dump(binned, f)
                    # del binned
                    # binned = {}

        update(bins, start_time, unsaved_docs, interval, i)
        bins.create_index([('bin', 'hashed')])
        bins.create_index([('name', 'hashed')])

        return bins


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
        for _bin in self.query_bins.find():
            __bin = _bin['bin']
            print('-------------------------------------------------')
            print('RESULTS FOR {}'.format(__bin))
            for result in self.binned.find({'bin':__bin}):
                print(result)
                names.append(result['name'])
                names.append(result['name'])

        print('Forward search done.')
        names = set(names)
        print(names)

        results = {}
        for name in names:
            print('-------------------------------------------------')
            print('Name: {}'.format(name))
            match = Match(name, self.query_bins, self.binned)
            match.find_edges()
            match.max_subgraph()
            match.plot_graph()
            # results[name] = []
            # print('searching {}'.format(name))
            # for _bin in self.binned.find({'name': name[0]}):
                # if _bin['idx1'] == name[1]:
                    # print('-------')
                    # print(_bin)
                    # for doc in self.query_bins.find({'bin':_bin['bin']}):
                        # print('MATCH:')
                        # results[name].append((doc['idx1'], doc['idx2']))
                        # print(doc)

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
            degrees=30, reset_querydb=True, dbname='test_bins')
    lookup.match()


def make_hash_table():
    print('Loading database and setting up lookup object...')
    # length cutoff of 2 turns or 10.8 angstroms
    lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            exposed_cutoff=0.3, length_cutoff=10.8, angstroms=5,
            degrees=30)
    print('Done.')
    # binned = lookup.bin_db(lookup.df)
    lookup.update_bin_db()
    # out = "binned_0p3/last.pkl"
    # with open(out, 'wb') as f:
        # pickle.dump(binned, f)

def make_test_hash_table():
    client = MongoClient()
    client['test_bins']['bins_5A_30D'].drop()
    lookup=HelixLookup(pd.read_pickle('out.pkl'), exposed_cutoff=0.3,
            length_cutoff=10.8, angstroms=5, degrees=30,
            dbname='test_bins')
    lookup.update_bin_db()


if __name__=='__main__':
    test()
    # make_hash_table()
    # make_test_hash_table()
