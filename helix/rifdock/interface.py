'''
Scan through PDB and find helices at protein interfaces where the
majority of the helix is part of the interface, with the idea that this
means the helix will be laying mostly flat along the protein, maximizing
its contact surface area.

Usage: 
    interface.py <num> [options]
    interface.py [options]

Options:
    --out=FOLDER, -o  Where to put the output  
    [default: interface_finder]
    --test=PDB  Test on a certain PDB
'''
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta import *
from pyrosetta.rosetta.utility import vector1_int
from rosetta.protocols.scoring import Interface
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import ScoreType
import pandas as pd
import docopt
import os, sys
import gzip


class InterfaceScore(object):
    def __init__(self, pose, dist=8.0):
        """'dist' is the maximum distance between two residues to
        consider them as part of an interface"""
        self.dist = dist
        self.pose = pose
        # ss_str = Dssp(self.pose).get_dssp_secstruct()
        # self.secstruct = contiguous_secstruct(ss_str)

        self.chains = 'A_BCDEFGH'

    def apply(self):

        # We are going to only change jump 1, so we know which jump to look at
        # in the Interface class
        movable_jumps = vector1_int(1)
        # Now we can alter the fold tree such that the jump is between the
        # chains we care about
        # chains = left_chains + '_' + right_chains
        sfxn = create_score_function('ref2015')
        sfxn(self.pose)
        print('SETTING UP FOLD TREE FOR INTERFACE {}'.format(self.chains))
        setup_foldtree(self.pose, self.chains, movable_jumps)


        # Set up interface class
        interface = Interface(1)
        interface.distance(self.dist)
        interface.calculate(self.pose)

        # Logic for filtering out residues not on chains of interest
        interface_list = []
        # interface_resis = []
        scoretypes = [ScoreType.fa_atr, ScoreType.fa_rep, 
                ScoreType.fa_sol, ScoreType.fa_elec, 
                ScoreType.hbond_bb_sc, ScoreType.hbond_sc]
        interface_score = 0
        for side in interface.pair_list():
            for resi in side:
                # Find out what chain the residue belongs to
                pdbinfo = self.pose.pdb_info().pose2pdb(resi)
                # resnum = pdbinfo.split(' ')[0]
                chain = pdbinfo.split(' ')[1]
                if chain == 'A':
                    resi_score = 0
                    for scoretype in scoretypes:
                        score = self.pose.energies().residue_total_energies(resi)[scoretype]
                        resi_score += score
                interface_score += resi_score
                # Find out what chain the resiude is interacting with
                closest = interface.closest_interface_residue(self.pose,
                        resi, self.dist)

        return interface_score


def test():
    init()
    idx = int(os.environ['SGE_TASK_ID']) - 1
    print('IDX = {}'.format(idx))
    num = int(args['<num>'])
    df = pd.DataFrame()
    start = idx * num
    stop = idx * num + num - 1
    print('START: {}'.format(start))
    print('STOP: {}'.format(stop))
    with gzip.open('nrpdb.gz', 'rb') as f:
        lines = f.readlines()[start:stop]

    errors = []
    df = pd.DataFrame()
    for line in lines:
        line = line.decode('utf-8')
        if not line.startswith('#'):
            try:
                print('Opening from line {}'.format(line))
                sys.stdout.flush()
                pdb = get_pdb_obj(str(line))
                if pdb:
                    pdb.interface_all_chains()
                    df = pd.concat([df, pdb.compile_helix_info()],
                            ignore_index=True)

            except Exception as e:
                print("Error scanning line: \n{}".format(line))
                print('Error was:')
                print(e)
                sys.stdout.flush()
                errors.append(line)
    outfolder = args['--out']
    df.to_csv('{}/pdb_interfaces_{}.csv'.format(outfolder, idx))
    df.to_pickle('{}/pdb_interfaces_{}.pkl'.format(outfolder, idx))

# Test of get_interface_residues(). Works, can delete.
def test(pdb):
    init()
    pose = pose_from_file(pdb)
    db = InterfaceScore(pose)
    print(db)
    print(db.apply())

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    if args['--test']:
        test(args['--test'])
    else:
        main()
