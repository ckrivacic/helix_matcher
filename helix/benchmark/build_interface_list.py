'''
Scan through PDB and find helices at protein interfaces where the
majority of the helix is part of the interface, with the idea that this
means the helix will be laying mostly flat along the protein, maximizing
its contact surface area.

Usage: 
    build_interface_list.py <num> [options]
    build_interface_list.py [options]

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
import pandas as pd
import docopt
import os, sys
import gzip


class PDB(object):
    def __init__(self, pdbid, dist=8.0):
        """'dist' is the maximum distance between two residues to
        consider them as part of an interface"""
        self.dist = dist
        self.pdbid = pdbid
        pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
        self.path = os.path.join(pdb_prefix, self.pdbid[1:3], 'pdb{}.ent.gz'.format(
            self.pdbid
            ))
        # self.path = pdbid
        self.pose = pose_from_file(self.path)
        ss_str = Dssp(self.pose).get_dssp_secstruct()
        self.secstruct = contiguous_secstruct(ss_str)

    def get_interface_residues(self, chains,
            name=None):

        # We are going to only change jump 1, so we know which jump to look at
        # in the Interface class
        movable_jumps = vector1_int(1)
        # Now we can alter the fold tree such that the jump is between the
        # chains we care about
        # chains = left_chains + '_' + right_chains
        sfxn = create_score_function('ref2015')
        sfxn(self.pose)
        weights = sfxn.weights()
        egraph = self.pose.energies().energy_graph()
        print('SETTING UP FOLD TREE FOR INTERFACE {}'.format(chains))
        setup_foldtree(self.pose, chains, movable_jumps)


        # Set up interface class
        interface = Interface(1)
        interface.distance(self.dist)
        interface.calculate(self.pose)

        if not name:
            name = self.pdbid

        # Logic for filtering out residues not on chains of interest
        interface_list = []
        contacts = interface.contact_list()
        # interface_resis = []
        for side in interface.pair_list():
            for resi in side:
                resi_score = 0
                # Find out what chain the residue belongs to
                pdbinfo = self.pose.pdb_info().pose2pdb(resi)
                # resnum = pdbinfo.split(' ')[0]
                chain = pdbinfo.split(' ')[1]
                # Find out what chain the resiude is interacting with
                closest = interface.closest_interface_residue(self.pose,
                        resi, self.dist)
                interacting_chain = self.pose.pdb_info().pose2pdb(closest).split(' ')[1]
                interacting_resi = self.pose.pdb_info().pose2pdb(closest).split(' ')[0]
                # If both this chain and its partner are in the chains we care
                # about, do something with them (in this case I'm spitting out a
                # string that will let me select all the interface residues in
                # PyMOL as a sanity check, you'll obviously want to save these
                # to a dictionary or something)
                if chain in list(chains) and interacting_chain in list(chains):
                    # interface_resis.append(resi)
                    for contact in contacts[resi]:
                        if self.pose.pdb_info().pose2pdb(contact).split(' ')[1] != chain:
                            edge = egraph.find_energy_edge(resi,
                                    contact).fill_energy_map()
                            filled = edge * weights
                            resi_score += filled.sum()
                    row = {'pdb': name,
                            'interface': chains,
                            'chain':chain,
                            'rosetta_resnum': resi,
                            'closest_chain': interacting_chain, 
                            'closest_rosetta_resnum': closest,
                            'residue_score': resi_score,
                            }
                            # 'closest_pdb_resnum': interacting_resi}
                            # 'pdb_resnum': resnum,
                    interface_list.append(row)

        return interface_list


    def interface_all_chains(self):
        pose_chains = []
        # Figure out what all the chains in the pose are
        for i in range(1, self.pose.num_chains() + 1):
            if self.pose.residue(self.pose.chain_begin(i)).is_protein():
                pdbinfo = self.pose.pdb_info().pose2pdb(self.pose.chain_begin(i))
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
            this_interface = self.get_interface_residues(interface_str)
            for row in this_interface:
                interfaces.append(row)

        # self.interfaces = pd.DataFrame(results)
        self.interfaces = pd.DataFrame(interfaces)
    
    def compile_helix_info(self):
        '''Make the helix dataframe'''
        rows = []
        for interface in set(self.interfaces['interface']):
            print(interface)
            this_interface = self.interfaces[self.interfaces['interface'] == interface]
            for helix in self.secstruct['H']:
                helix_interface_resis = []
                helix_resis = [x for x in range(helix[0], helix[1] + 1)]
                for n in helix_resis:
                    if n in this_interface['rosetta_resnum'].tolist():
                        helix_interface_resis.append(n)
                interacting_chains = []
                interface_score = 0
                interface_scores = []
                interacting_residues_by_score = 0
                for resi in helix_interface_resis:
                    this_residue = this_interface[this_interface['rosetta_resnum'] == resi]
                    interacting_chains.append(this_residue['closest_chain'].tolist()[0])
                    interface_score += this_residue['residue_score'].tolist()[0]
                    interface_scores.append(this_residue['residue_score'].tolist()[0])
                    if this_residue['residue_score'].tolist()[0] < -1.5:
                        interacting_residues_by_score += 1

                chain = self.pose.pdb_info().pose2pdb(helix[0]).split(' ')[1]
                pdb_resis = [self.pose.pdb_info().pose2pdb(x).split(' ')[0] for x in
                        helix_resis]

                if len(interacting_chains) > 0:
                    interacting_chain = max(set(interacting_chains),
                            key=interacting_chains.count)
                else:
                    interacting_chain = None
                row = {
                        'name': self.pdbid,
                        'interface': interface,
                        'chain': chain,
                        'rosetta_resis': helix_resis,
                        'interface_resis': helix_interface_resis,
                        'pdb_resis': pdb_resis,
                        'interacting_chain': interacting_chain,
                        'all_interacting_chains':
                        interacting_chains,
                        'interacting_length':
                        len(helix_interface_resis),
                        'interface_score': interface_score,
                        'interacting_residues_by_score':
                        interacting_residues_by_score,
                        'interface_scores': interface_scores,
                        }

                rows.append(row)

        return pd.DataFrame(rows)



def contiguous_secstruct(ss_str):
	"""
    This function takes the output from rosetta.core.scoring.dssp.Dssp(pair.pose), which is
    a string of H, L, and E's that denote secondary struct at each residue.

    Returns a dictionary. Keys = H, L, or E. Values = lists of tuples, each corresp to one continuous secstruct element
    """
    ss_positions = {}
    ss_positions['H'] = []
    ss_positions['L'] = []
    ss_positions['E'] = []

	start = 0
	for i in range(0, len(ss_str)):
		if ss_str[i] != ss_str[start]:
			ss_positions[ss_str[start]].append((start + 1, i))
			start = i

	if i + 1== len(ss_str):
		ss_positions[ss_str[start]].append((start + 1, i + 1))

	return ss_positions


def get_pdb_obj(line):
    pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
    fields = list(filter(None, line.split(' ')))
    pdb = fields[0].lower()
    chain = fields[1]
    rep = int(fields[5])
    if rep == 1:
        return PDB(pdb)
    else:
        return None


def main():
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
    db = PDB(pdb.lower())
    db.interface_all_chains()
    df = db.compile_helix_info()
    print(df)

if __name__=='__main__':
    args = docopt.docopt(__doc__)
    if args['--test']:
        test(args['--test'])
    else:
        main()
