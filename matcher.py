import pandas as pd
import numpy as np
from roseasy.utils import numeric
from itertools import product
import os, psutil, sys
import pickle
import subprocess
from pymongo import MongoClient
from pyrosetta import init, pose_from_file
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


def bin_array(array, bins):
    '''Digitize a numpy array'''
    inds = np.digitize(array, bins)
    # print(inds)
    binned = tuple([bins[inds[n]-1] for n in range(array.size)])
    # print(binned)
    return binned


class HelixLookup(object):

    def __init__(self, df, exposed_cutoff=0.5, length_cutoff=10.8,
            query_df=None, query_name=None, angstroms=2, degrees=20):
        # Setup pymongo
        # self.client = MongoClient()
        self.client = MongoClient()
        try:
            print(self.client.server_info())
        except:
            print('Problem with MongoDB server connection; attempting to open...')
            subprocess.Popen(['mongod', '--dbpath',
                '/Users/codykrivacic/data/db'])
            print(self.client.server_info())


        self.df = df
        self.df['idx'] = self.df.index
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
        self.binned = self.client['helix_bins'][binned_name]

        if query_df is not None:
            self.query_df = query_df
            if query_name is None:
                query_name = self.query_df.iloc[0]['name']
            self.query_df['idx'] = self.query_df.index
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
        self.bin_db(self.df, 'helix_bins', check_name=True)

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
                    idx1 = combination[0]['idx']
                    idx2 = combination[1]['idx']
                    transform = numeric.Transformation(combination[0]['vector'],
                            combination[1]['vector'])
                    rot = R.from_matrix(transform.rotation)
                    rot = rot.as_euler('xyz', degrees=True)
                    # combination[0]['rotation'] = rot
                    # combination[0]['translation'] = transform.translation

                    rbin = bin_array(rot, self.rbins)
                    tbin = bin_array(transform.translation, self.tbins)
                    rbin2 = bin_array(rot, self.rbins + (self.degrees/2))
                    tbin2 = bin_array(transform.translation, self.tbins +
                            (self.angstroms/2))

                    x = [tbin[0], tbin2[0]]
                    y = [tbin[1], tbin2[1]]
                    z = [tbin[2], tbin2[2]]
                    phi = [rbin[0], rbin2[0]]
                    psi = [rbin[1], rbin2[1]]
                    om = [rbin[2], rbin2[2]]

                    for bin_12 in product(x, y, z, phi, psi,
                        om):
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
            print('DICT FOR {}'.format(__bin))
            for result in self.binned.find({'bin':__bin}):
                names.append(result['name'])

        results = {}
        for name in names:
            results[name] = []
            for _bin in self.binned.find({'name': name}):
                for doc in self.query_bins.find({'bin':_bin['bin']}):
                    results[name].append(self.query_bins['bin'])

        for key in results:
            print('PDB {} had {} matching results'.format(
                key, len(results[key])
                ))


def test():
    import scan_helices

    test_path = 'test_files/6r9d.cif'
    init()
    pose = pose_from_file(test_path).split_by_chain(1)
    print(pose.size())
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices()
    helices = pd.DataFrame(helices)
    print(helices)
    helices = helices[helices['percent_exposed'] > 0.3]
    print(helices)

    # lookup = HelixLookup(pd.read_pickle('dataframes/final.pkl'),
            # query_df=helices, query_name='6r9d')
    lookup = HelixLookup(pd.DataFrame(),
            query_df=helices, query_name='6r9d', angstroms=5, degrees=30)
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


if __name__=='__main__':
    test()
    # make_hash_table()
