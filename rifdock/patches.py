from pyrosetta.rosetta.core.select import residue_selector
import sys, os
import pandas as pd
import numpy as np
from pyrosetta import *
sys.path.insert(1,
        os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')))
from numeric import *
from utils import *

class Patches(object):
    '''
    Class which holds information about protein interface patches
    '''
    def __init__(self, pose):
        if pose.num_chains() > 1:
            print("WARNING: Pose has more than one chain. Interface "\
                    "residues between chains may not be counted as exposed,"\
                    " and therefore won't be used to determine surface "\
                    "patches. For better results, pass a single chain "\
                    "using pose.split_by_chain(chain_num)")
        self.pose = pose
        self.reslist = None
        self.pdb_reslist = None

    def set_reslist(self, reslist):
        '''
        Give the object your own reslist. Useful if you want to only
        align to interfaces.
        Make sure to run determine_surface_residues if you want to
        narrow down this reslist to solvent-exposed residues.
        '''
        if type(reslist[0])==type(True):
            self.reslist = res_selector_to_size_list(reslist)
        else:
            self.reslist = intlist_to_vector1_size(reslist) 
        # Use PDB numbering so that we can split chains and still
        # get the correct residues.
        self.pdb_reslist = reslist_to_pdb_numbers(self.reslist,
                self.pose)

    def determine_surface_residues(self):
        """Restrict reslist to surface residues"""
        # Use only the chain of the given reslist if possible
        if self.reslist:
            chain = self.pose.chain(self.reslist.front())
            chain_pose = self.pose.split_by_chain(chain)
        else:
            chain_pose = self.pose


        surface_selector = residue_selector.LayerSelector()
        surface_selector.set_layers(False, False, True)
        #try:
        print(surface_selector.apply(chain_pose))
        print(res_selector_to_size_list(surface_selector.apply(chain_pose)))
        reslist =\
                res_selector_to_size_list(surface_selector.apply(chain_pose))
        reslist = correct_resnums(chain_pose, reslist, self.pose)
        #except:
        #    reslist = []
        if self.reslist:
            # If a reslist is already defined, only  take surface
            # residues that are in that reslist
            reslist = [res for res in reslist if res in self.reslist]
        self.reslist = reslist

    def map_residues(self):
        '''
        Makes a dictionary that stores pairwise distances between all
        residues in reslist
        '''
        resmap = {}
        for res1 in self.reslist:
            if res1 not in resmap:
                resmap[res1] = {}
            for res2 in self.reslist:
                # Don't calculate twice
                if res2 in resmap and res1 in resmap[res2]: 
                    resmap[res1][res2] = resmap[res2][res1]
                else:
                    xyz1 = self.pose.residue(res1).xyz('CA')
                    xyz2 = self.pose.residue(res2).xyz('CA')
                    resmap[res1][res2] = euclidean_distance(xyz1, xyz2)
        if len(resmap) > 0:
            resmap = pd.DataFrame(resmap).fillna(0).unstack().reset_index()
            resmap.columns = ['res1', 'res2', 'dist']
            self.resmap = resmap
        else:
            self.resmap = None

    def nearest_n_residues(self, resnum, n, cutoff=30.0, pymol=False):
        if np.any(self.resmap) and not self.resmap.empty:
            neighbors = self.resmap[(self.resmap['res1']==resnum) &
                    (self.resmap['dist'] <
                        cutoff)].sort_values(by='dist')['res2']
            if not pymol:
                return set(neighbors[0:n].tolist())
            else:
                return reslist_to_pdb_numbers(set(neighbors[0:n].tolist()))
        else:
            return None


if __name__=='__main__':
    init()
    pdbid = sys.argv[1]
    pose = pose_from_rcsb(pdbid, 'test_inputs')

    patches = Patches(pose)
    #patches.set_reslist([7,11,14, 8])
    patches.determine_surface_residues()
    print('reslist:')
    print(patches.reslist)
    
    patches.map_residues()
    print(patches.resmap)
    #print(patches.resmap)
    print(patches.nearest_n_residues(11, 10))
