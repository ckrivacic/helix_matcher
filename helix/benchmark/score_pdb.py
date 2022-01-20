'''
Find average score of a number of structures in the PDB for each buried
amino acid. Possibly all interface amino acids as well.

Usage:
    score_pdb.py <nstruct> [options]

Options:
    --test=PDB  Test on a certain pdb
'''
import docopt
from pyrosetta import *
from pyrosetta.rosetta.core.scoring import ScoreType
from helix.utils import utils
import pandas as pd
import os


class PDBInterface(object):
    def __init__(pdbpath, sfxn=None):
        if not sfxn:
            self.sfxn = create_score_function('ref2015')
        else:
            self.sfxn = sfxn
        self.pose = pose_from_file(pdbpath)
        self.ss_str = Dssp(self.pose).get_dssp_secstruct()

    def get_scores(interface, pdbid):
        '''For a residue, get the cross-chain scores'''

        # We are going to only change jump 1, so we know which jump to look at
        # in the Interface class
        movable_jumps = vector1_int(1)
        # Now we can alter the fold tree such that the jump is between the
        # chains we care about
        # chains = left_chains + '_' + right_chains
        self.sfxn(self.pose)
        weights = self.sfxn.weights()
        egraph = self.pose.energies().energy_graph()
        print('SETTING UP FOLD TREE FOR INTERFACE {}'.format(chains))
        setup_foldtree(self.pose, chains, movable_jumps)

        scoretypes = [ScoreType.fa_atr, ScoreType.fa_rep, 
                ScoreType.fa_sol, ScoreType.fa_elec, 
                ScoreType.hbond_bb_sc, ScoreType.hbond_sc]

        egraph = self.pose.energies().energy_graph()
        contacts = interface.contact_list()
        interface_list = []
        # for side in interface.pair_list():
        side = interface.pair_list()[0]
        for resi in side:
            pdbinfo = self.pose.pdb_info().pose2pdb(resi)
            chain = pdbinfo.split(' ')[1]
            row = {
                    'resnum': resi,
                    'chain': chain,
                    'pdb': pdbid,
                    'restype': self.pose.residue(resi).name3()
                    }
            # Fill scoretype columns
            for scoretype in scoretypes:
                row[str(st).split('.')[1]] = 0
            # Iterate through contacts, adding scoretypes and getting
            # total cross-chain energy
            for contact in contacts[resi]:
                edge = egraph.find_energy_edge(resi,
                        contact).fill_energy_map()
                filled = edge * weights
                for scoretype in scoretypes:
                    sc = edge.get(scoretype)
                    row[str(scoretype).split('.')[1]] += sc
            row['total_crosschain'] = filled.sum()
            row['secstruct'] = self.ss_str[resi - 1]
            interface_list.append(row)

        return interface_list


    def interface_all_chains(self):
        pose_chains = []
        # Figure out what all the chains of the pose are
        for i in range(1, self.pose.num_chains() + 1):
            if self.pose.residue(self.pose.chain_begin(i)).is_protein():
                pdbinfo =
                self.pose.pdb_info().pose2pdb(self.pose.chain_begin(i))
                chain = pdbinfo.split(' ')[1]
                pose_chains.append(chain)


        interfaces = []
        pose_chains = list(set(pose_chains))
        # Now get interface residues using each chain as the query
        for i in range(0, len(pose_chains)):
            query_chain = pose_chains[i]
            left = pose_chains[0:i]
            right = pose_chains[i+1:len(pose_chains)]
            # Left and right are used differently here than in the
            # get_interface function; here it's just the chains to the left
            # and right of the query chain in this list. We combine them,
            # and that whole thing becomes the right in the get_interface
            # function.
            other_chains = ''.join(left + right)
            interface_str = '{}_{}'.format(query_chain, other_chains)
            this_interface = self.gget_scores(interface_str)
            for row in this_interface:
                interfaces.append(row)

        self.interfaces = pd.DataFrame(interfaces)


def main():
    init()
    args = docopt.docopt(__doc__)
    if args['--test']:
        pdbid = args['--test']
        pdbpath = os.path.join('/wynton/home/database/pdb/remediated/pdb', 
                pdbid[1:3], 'pdb{}.ent.gz'.format(pdbid))
        pdbpaths = [pdbpath]
    else:
        idx = int(os.environ['SGE_TASK_ID']) - 1
        num = int(args['<num>'])
        start = idx * num
        stop = idx * num -+ num - 1
        
        all_pdbpaths = utils.get_all_pdbs(pdbid=False)
        if stop > len(all_pdbids):
            stop = len(all_pdbids)
        pdbpaths = all_pdbpaths[start:stop]

    df = pd.DataFrame()

    errors = []
    sfxn = create_score_function('ref2015')
    for pdbpath in pdbpaths:
        try:
            pdb_obj = PDBInterface(pdbpath, sfxn=sfxn)
            df = pd.concat([df, pdb_obj.interface_all_chains()],
                    ignore_index=True)
        except Exception as e:
            print('Error analyzing interface of {}'.format(pdbpath))
            print('Error was:')
            print(e)

    outfolder = 'residue_scores/'
    df.to_pickle(os.path.join(outfolder,
        'pdb_interface_scores_{}.pkl'.format(idx)))
