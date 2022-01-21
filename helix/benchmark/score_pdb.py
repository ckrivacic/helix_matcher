'''
Find score (and each scoreterm) for all interface residues in the PDB.

Usage:
    score_pdb.py <nstruct> [options]

Options:
    --test=PDB  Test on a certain pdb
    --local, -l  Run a test run off wynton
'''
import docopt
import traceback
from pyrosetta import *
from pyrosetta.rosetta.utility import vector1_int
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.protocols.scoring import Interface
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from helix.utils import utils
import pandas as pd
import os


class PDBInterface(object):
    def __init__(self, pdbpath, sfxn=None, pdbid=None, dist=8.0):
        if not sfxn:
            self.sfxn = create_score_function('ref2015')
        else:
            self.sfxn = sfxn
        self.pose = pose_from_file(pdbpath)
        minmover_str = '''
        <MinMover name='minimize' jump='all' bb='true' chi='true'/>
        '''
        minmover = XmlObjects.static_get_mover(minmover_str)
        minmover.apply(self.pose)
        self.ss_str = Dssp(self.pose).get_dssp_secstruct()
        self.dist = dist
        if pdbid:
            self.pdbid = pdbid
        else:
            self.pdbid='0000'

        self.setup_selectors()

    def setup_selectors(self):
        '''Setup some layer selectors and determine boundary and buried
        residues'''
        buried = residue_selector.LayerSelector()
        buried.set_layers(True, False, False)
        boundary = residue_selector.LayerSelector()
        boundary.set_layers(False, True, False)

        buried_selector = buried.apply(self.pose)
        boundary_selector = boundary.apply(self.pose)

        self.buried = utils.res_selector_to_size_list(buried_selector,
                pylist=True)
        self.boundary = utils.res_selector_to_size_list(boundary_selector, 
                pylist=True)

    def get_scores(self, chains):
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
        tries = 0
        while tries < 3:
            # It's possible to have the jump go to a ligand, which would
            # cause the Interface to look for ligand-interacting
            # residues.
            # Can be avoided by choosing a random residue to jump to for
            # the 2nd partner (rather than COM), but it could still go to a ligand,
            # so I'll give it a few tries to be safe. Yes I know how
            # ugly this is.
            setup_foldtree(self.pose, chains, movable_jumps, True)
            # print('NEW FOLD TREE')
            # print(self.pose.fold_tree())

            scoretypes = [ScoreType.fa_atr, ScoreType.fa_rep, 
                    ScoreType.fa_sol, ScoreType.fa_elec, 
                    ScoreType.hbond_bb_sc, ScoreType.hbond_sc,
                    ]
            scoretypes_full = scoretypes.copy()
            scoretypes_full.extend([
                ScoreType.pro_close, ScoreType.fa_intra_rep,
                ScoreType.dslf_fa13, ScoreType.rama, ScoreType.omega,
                ScoreType.fa_dun, ScoreType.p_aa_pp, ScoreType.ref,
                ScoreType.hbond_sr_bb
                ])

            egraph = self.pose.energies().energy_graph()

            # Set up interface class
            interface = Interface(1)
            interface.distance(self.dist)
            interface.calculate(self.pose)
            contacts = interface.contact_list()
            interface_list = []

            # Get all residues in interface from the chain of interest
            pdbinf = self.pose.pdb_info()
            side = []
            for intside in interface.pair_list():
                for resi in intside:
                    if pdbinf.pose2pdb(resi).split(' ')[1] ==\
                            chains.split('_')[0] and\
                            self.pose.residue(resi).is_protein():
                        side.append(resi)
            print('SIDE')
            print(side)
            if len(side) > 0:
                break
            tries += 1

        # side = interface.pair_list()[0]
        for resi in side:
            pdbinfo = self.pose.pdb_info().pose2pdb(resi)
            chain = pdbinfo.split(' ')[1]
            row = {
                    'resnum': resi,
                    'pdb_resnum': pdbinfo.split(' ')[0],
                    'chain': chain,
                    'pdb': self.pdbid,
                    'restype': self.pose.residue(resi).name3()
                    }
            # Fill scoretype columns
            for scoretype in scoretypes:
                row[str(scoretype).split('.')[1] + '_cc'] = 0
            # Iterate through contacts, adding scoretypes and getting
            # total cross-chain energy
            total_crosschain = 0
            for contact in contacts[resi]:
                edge = egraph.find_energy_edge(resi,
                        contact).fill_energy_map()
                filled = edge * weights
                total_crosschain += filled.sum()
                for scoretype in scoretypes:
                    sc = edge.get(scoretype)
                    row[str(scoretype).split('.')[1] + '_cc'] += sc
            row['total_crosschain'] = total_crosschain
            for scoretype in scoretypes_full:
                sc = self.pose.energies().residue_total_energies(resi).get(scoretype)
                row[str(scoretype).split('.')[1] + '_tot'] = sc
            row['total_energy'] = self.pose.energies().residue_total_energy(resi)
            if resi in self.buried:
                row['burial'] = 'buried'
            elif resi in self.boundary:
                row['burial'] = 'boundary'
            else:
                row['burial'] = 'surface'
            row['secstruct'] = self.ss_str[resi - 1]
            if len(contacts[resi]) > 0:
                row['contacts'] = len(contacts[resi])
                interface_list.append(row)

        return interface_list


    def interface_all_chains(self):
        pose_chains = []
        # Figure out what all the chains of the pose are
        for i in range(1, self.pose.num_chains() + 1):
            if self.pose.residue(self.pose.chain_begin(i)).is_protein():
                pdbinfo = self.pose.pdb_info().pose2pdb(self.pose.chain_begin(i))
                chain = pdbinfo.split(' ')[1]
                pose_chains.append(chain)


        interfaces = []
        pose_chains = list(set(pose_chains))
        if len(pose_chains) < 2:
            self.interfaces = pd.DataFrame()
            return
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
            this_interface = self.get_scores(interface_str)
            for row in this_interface:
                interfaces.append(row)

        self.interfaces = pd.DataFrame(interfaces)
        return self.interfaces


def main():
    init('-ignore_unrecognized_res')
    args = docopt.docopt(__doc__)
    if args['--test']:
        pdbid = args['--test']
        if not args['--local']:
            pdbpath = os.path.join('/wynton/home/database/pdb/remediated/pdb', 
                    pdbid[1:3], 'pdb{}.ent.gz'.format(pdbid))
        else: 
            import wget
            url = 'http://files.rcsb.org/download/' + pdbid + '.pdb'
            if not os.path.exists(pdbid + '.pdb'):
                wget.download(url, pdbid + '.pdb')
            pdbpath = pdbid + '.pdb'
        idx = 1
        pdbpaths = [pdbpath]
    else:
        idx = int(os.environ['SGE_TASK_ID']) - 1
        num = int(args['<nstruct>'])
        start = idx * num
        stop = idx * num + num - 1
        
        all_pdbpaths = utils.get_pdb_list(pdbid=False)
        if stop > len(all_pdbpaths):
            stop = len(all_pdbpaths)
        pdbpaths = all_pdbpaths[start:stop]

    df = pd.DataFrame()

    errors = []
    sfxn = create_score_function('ref2015')
    for pdbpath in pdbpaths:
        if not args['--test']:
            pdbid = os.path.basename(pdbpath)[3:7]
        try:
            pdb_obj = PDBInterface(pdbpath, sfxn=sfxn)
            df = pd.concat([df, pdb_obj.interface_all_chains()],
                    ignore_index=True)
        except Exception as e:
            print('Error analyzing interface of {}'.format(pdbpath))
            print('Error was:')
            print(e)
            print(traceback.format_exc())

    outfolder = 'residue_scores_min/'
    if args['--test']:
        df.to_csv(os.path.join(outfolder,
            'min_pdb_interface_scores_{}.csv'.format(idx)))
        df.to_pickle(os.path.join(outfolder,
            'min_pdb_interface_scores_{}.pkl'.format(idx)))
    else:
        df.to_pickle(os.path.join(outfolder,
            'min_pdb_interface_scores_{}.pkl'.format(idx)))

if __name__=='__main__':
    main()
