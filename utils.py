from pyrosetta import *
from pyrosetta.rosetta.core.select import residue_selector
from numeric import intlist_to_vector1_size
import os, wget

def three_to_one(restype):
    three_to_one = {
            'ala':'a', 'arg':'r', 'asn':'n', 'asp':'d', 'cys':'c',
            'glu':'e', 'gln':'q', 'gly':'g', 'his':'h', 'ile':'i',
            'leu':'l', 'lys':'k', 'met':'m', 'phe':'f', 'pro':'p',
            'ser':'s', 'thr':'t', 'trp':'w', 'tyr':'y', 'val':'v'
            }

    return three_to_one[restype.lower()].upper()

def pose_from_rcsb(pdbid, prefix=None):
    if prefix:
        path = os.path.join(prefix,pdbid)
    else:
        path = pdbid
    if not os.path.exists(path + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    pyrosetta.toolbox.cleanATOM(path + '.pdb')
    pose = rosetta.core.import_pose.get_pdb_and_cleanup(path + '.clean.pdb')

    return pose


def list_to_str(l):
    return ','.join(list(str(i) for i in l))


def list_to_res_selector(l):
    return residue_selector.ResidueIndexSelector(list_to_str(l))


def res_selector_to_size_list(resselector):
    size_list = []
    for i, boolean in enumerate(resselector):
        if boolean == True:
            size_list.append(int(i + 1))

    return intlist_to_vector1_size(size_list)


def reslist_to_pdb_numbers(reslist, pose, chain=None):
    poselist = []
    for res in reslist:
        poselist.append(pose.pdb_info().pose2pdb(res))

    return poselist


def int_list_to_pdb_numbers(reslist, chain='Z'):
    '''
    Takes a list of residues with PDB numbering and formats it as a list
    of strings with the chain number included (i.e., '19_Z')
    '''
    pdblist = []
    for res in reslist:
        if type(res) == type('hi'):
            spl = res.split(' ')
            resnum = res.split(' ')[0]
            pdblist.append(' '.join([str(resnum), chain]))
        else:
            pdblist.append(' '.join([str(res), chain]))

    return pdblist


def rosetta_numbers_from_pdb(reslist, pose, chain='A'):
    '''Get list of rosetta numbers from a list of PDB numbers'''
    rosetta_list = []
    for resi in reslist:
        if type(resi) == type('hi'):
            spl = resi.split(' ')
            rosetta_list.append(pose.pdb_info().pdb2pose(spl[1], spl[0]))
        else:
            rosetta_list.append(pose.pdb_info().pdb2pose(chain, resi))

    return rosetta_list


def reslist_to_selstr(reslist, chain='A'):
    selection = []
    for res in reslist:
        try:
            splt = res.split(' ')
            resnum = splt[0]
            chain = splt[1]
        except:
            # Default to chain A if given a list of integers instead of
            # resnum and chain. Can add an option to specify chain
            # later.
            resnum = res
            chain = chain
        selstr = ' (resi {} and chain {}) '.format(resnum, chain)
        selection.append(selstr)

    return 'or'.join(selection)


def reslist_to_selstr_chain(reslist, chain=None):
    selection = []
    for res in reslist:
        if not chain:
            splt = res.split(' ')
            resnum = splt[0]
            chain = splt[1]
        else:
            # Default to chain A if given a list of integers instead of
            # resnum and chain. Can add an option to specify chain
            # later.
            resnum = res.split(' ')[0]
            chain = chain
        selstr = ' (resi {} and chain {}) '.format(resnum, chain)
        selection.append(selstr)

    return 'or'.join(selection)


def correct_resnums(initial_pose, reslist, final_pose):
    """Get rosetta numbering for final_pose from a reslist for the
    intial_pose"""
    corrected_residues = []
    for res in reslist:
        pdbnum = initial_pose.pdb_info().pose2pdb(res)
        pdbres = int(pdbnum.split(' ')[0])
        pdbchain = pdbnum.split(' ')[1]
        rosettanum = final_pose.pdb_info().pdb2pose(pdbchain, pdbres)
        corrected_residues.append(rosettanum)

    return corrected_residues

# Write a test for correct_resnums
def test_correct_resnums():
    pose = pose_from_file('test_inputs/3r7c.clean.pdb')
    chain_pose = pose.split_by_chain(2)
    reslist = [60]
    correct_resnums(pose, reslist, pose)

if __name__=='__main__':
    init()
    test_correct_resnums()
