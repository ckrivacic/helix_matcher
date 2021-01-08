import pandas as pd
import numpy as np
from roseasy.utils import numeric
from itertools import product
import os, psutil
import pickle
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

    def __init__(self, path, exposed_cutoff=0.5, binned_path=None):
        self.db_path = path
        self.df = pd.read_pickle(path)
        self.df['idx'] = self.df.index
        self.df = self.df[self.df['percent_exposed'] > exposed_cutoff]
        self.degrees = 20
        self.angstroms = 2
        self.setup_bins()
        if binned_path:
            with open(binned_path, 'rb') as f:
                self.binned = pickle.load(f)

    def setup_bins(self):
        nrbins = int(360//self.degrees) + 1
        self.rbins = np.linspace(-180, 180, nrbins)
        tstart = -10000
        tstop = 10000
        ntbins = int((tstop - tstart) // self.angstroms) + 1
        self.tbins = np.linspace(tstart, tstop, ntbins)

    def bin_db(self):
        '''
        Bins only need to be the lengths of the two ends of the helices
        from one another...
        '''

        from scipy.spatial.transform import Rotation as R
        bin_size = 1
        self.binned = {}
        i = 0
        j = 0
        for name, group in self.df.groupby(['name']):
            i += 1
            if i%1000 == 0:
                mem_used = psutil.Process(os.getpid()).memory_info().rss
                print('{} PDBs processed so far.'.format(i))
                print('Currently using {} G of memory'.format(
                    mem_used * 10**-9
                    ))
                # If more than 6G memory used, dump just in case.
                if mem_used > 6 * 10**9:
                    out = "binned/{}.pkl".format(j)
                    j += 1
                    with open(out, 'wb') as f:
                        pickle.dump(self.binned, f)
                    self.binned = {}

            for combination in product(group.T.to_dict().values(),
                    repeat=2):
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
                if (tbin, rbin) not in self.binned:
                    self.binned[(tbin, rbin)] = {}
                if name not in self.binned[(tbin, rbin)]:
                    self.binned[(tbin, rbin)][name] = []
                self.binned[(tbin,
                    rbin)][name].append((idx1, idx2))

                rbin = bin_array(rot, self.rbins + (self.degrees/2))
                tbin = bin_array(transform.translation, self.tbins +
                        (self.angstroms/2))

                if (tbin, rbin) not in self.binned:
                    self.binned[(tbin, rbin)] = {}
                if name not in self.binned[(tbin, rbin)]:
                    self.binned[(tbin, rbin)][name] = []
                self.binned[(tbin,
                    rbin)][name].append((idx1, idx2))



def test():
    lookup = HelixLookup('dataframes/final.pkl',
            binned_path='binned/final.pkl')
    # lookup.bin_db()
    # out = "binned/last.pkl"
    # with open(out, 'wb') as f:
        # pickle.dump(lookup.binned, f)


if __name__=='__main__':
    test()
