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
        weights = sfxn.weights()
        sfxn(self.pose)
        print('SETTING UP FOLD TREE FOR INTERFACE {}'.format(self.chains))
        setup_foldtree(self.pose, self.chains, movable_jumps)


        # Set up interface class
        interface = Interface(1)
        interface.distance(self.dist)
        interface.calculate(self.pose)

        egraph = self.pose.energies().energy_graph()
        hbondset = bindings.pose.get_hbonds(self.pose)

        # Logic for filtering out residues not on chains of interest
        # interface_list = []
        contacts = interface.contact_list()
        # interface_resis = []
        scoretypes = [ScoreType.fa_atr, ScoreType.fa_rep, 
                ScoreType.fa_sol, ScoreType.fa_elec, 
                ScoreType.hbond_bb_sc, ScoreType.hbond_sc]
        interface_score = 0
        self.n_hbonds = 0
        for side in interface.pair_list():
            for resi in side:
                resi_score = 0
                # Find out what chain the residue belongs to
                pdbinfo = self.pose.pdb_info().pose2pdb(resi)
                # resnum = pdbinfo.split(' ')[0]
                chain = pdbinfo.split(' ')[1]
                if chain == 'B':
                    for contact in contacts[resi]:
                        if self.pose.pdb_info().pose2pdb(contact).split(' ')[1] != chain:
                            edge = egraph.find_energy_edge(resi,
                                    contact).fill_energy_map()
                            filled = edge * weights
                            resi_score += filled.sum()
                    for hbond in hbondset.residue_hbonds(resi):
                        if hbond.energy() < -0.5 and not\
                                is_same_chain(self.pose, resi, hbond):
                            self.n_hbonds += 1

                    # for scoretype in scoretypes:
                        # score = self.pose.energies().residue_total_energies(resi)[scoretype]
                        # if scoretype ==\
                                # ScoreType.hbond_bb_sc or scoretype\
                                # == ScoreType.hbond_sc:
                            # if score < -0.5:
                                # self.n_hbonds += 1
                        # resi_score += score
                interface_score += resi_score
                # Find out what chain the resiude is interacting with
                # closest = interface.closest_interface_residue(self.pose,
                        # resi, self.dist)

        return interface_score


def is_same_chain(pose, resi, hbond):
    chain1 = pose.pdb_info(resi).pose2pdb(hbond.acc_res).split(' ')[1]
    chain2 = pose.pdb_info(resi).pose2pdb(hbond.don_res).split(' ')[1]
    return chain1 == chain2


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


if __name__=='__main__':
    args = docopt.docopt(__doc__)
    if args['--test']:
        test(args['--test'])
    else:
        main()
