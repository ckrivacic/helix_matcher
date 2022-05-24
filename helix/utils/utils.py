from pyrosetta import *
from pyrosetta.rosetta.core.select import residue_selector
from helix.utils.numeric import intlist_to_vector1_size
import os, wget, sys
import networkx as nx
import prody


def safe_load(pickledf):
    import pandas as pd
    try:
        dataframe = pd.read_pickle(pickledf)
    except:
        import pickle5 as pickle
        with open(pickledf, 'rb') as f:
            dataframe = pickle.load(f)

    return dataframe


def get_pdb_list(pdbid=True):
    import glob
    prefix = '/wynton/home/database/pdb/remediated/pdb/'
    all_pdbs = glob.glob(os.path.join(prefix, '*', '*.ent.gz'))
    if pdbid:
        pdbids = []
        for f in all_pdbs:
            pdbids.append(os.path.basename(f)[3:7])
        return sorted(pdbids)
    else:
        return sorted(all_pdbs)


def three_to_one(restype):
    three_to_one = {
            'ala':'a', 'arg':'r', 'asn':'n', 'asp':'d', 'cys':'c',
            'glu':'e', 'gln':'q', 'gly':'g', 'his':'h', 'ile':'i',
            'leu':'l', 'lys':'k', 'met':'m', 'phe':'f', 'pro':'p',
            'ser':'s', 'thr':'t', 'trp':'w', 'tyr':'y', 'val':'v'
            }

    return three_to_one[restype.lower()].upper()

def pose_from_rcsb(pdbid, prefix=None, use_prody=False):
    if prefix:
        path = os.path.join(prefix,pdbid)
    else:
        path = pdbid
    if not os.path.exists(path + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')
    if use_prody:
        return prody.parsePDB(path + '.pdb')
    else:
        pyrosetta.toolbox.cleanATOM(path + '.pdb')
        pose = rosetta.core.import_pose.get_pdb_and_cleanup(path + '.clean.pdb')
        return pose


def pose_from_wynton(pdbid,
        prefix='/wynton/home/database/pdb/remediated/pdb/', clean=False,
        use_prody=False, header=False):
    '''
    Import a pose from the Wynton database. Assumes PDB is stored with
    the format <base_pdb_folder>/<pdbid[1:3]>/pdb<pdbid>.ent.gz
    '''
    from pyrosetta import toolbox
    pdbid = pdbid.lower()
    folder = os.path.join(prefix, pdbid[1:3])
    fname = 'pdb{}.ent.gz'.format(pdbid)
    fpath = os.path.join(folder, fname)

    print('HELIX opening following from Wynton:')
    print(pdbid)
    print(fpath)

    if clean:
        print('Cleaning PDB file')
        toolbox.cleanATOM(fpath, out_file=f'{pdbid}.clean.pdb')
        return  pose_from_file(f'{pdbid}.clean.pdb')
    elif use_prody:
        return prody.parsePDB(fpath, header=header)
    else:
        return pose_from_file(fpath)


def download_and_clean_pdb(pdbid, prefix=None):
    if prefix:
        path = os.path.join(prefix, pdbid)
    else:
        path = pdbid
    if not os.path.exists(path + '.pdb'):
        url = 'https://files.rcsb.org/download/' + pdbid + '.pdb'
        wget.download(url, path + '.pdb')

    pyrosetta.toolbox.cleanATOM(path + '.pdb')
    os.remove(path + '.pdb')

    return path + '.clean.pdb'


def list_to_str(l):
    return ','.join(list(str(i) for i in l))


def list_to_res_selector(l):
    return residue_selector.ResidueIndexSelector(list_to_str(l))


def res_selector_to_size_list(resselector, pylist=False):
    size_list = []
    for i, boolean in enumerate(resselector):
        if boolean == True:
            size_list.append(int(i + 1))

    if pylist:
        return size_list
    else:
        return intlist_to_vector1_size(size_list)


def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector


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


def max_subgraph(graph):
    max_subgraph_len = 0
    max_subgraph = None
    subgraphs = []
    for f in nx.find_cliques(graph):
        if len(f) > max_subgraph_len:
            subgraphs = []
            max_subgraph_len = len(f)
            subgraphs.append(f)
        elif len(f) == max_subgraph_len:
            subgraphs.append(f)

    return subgraphs


def run_command(cmd, environment=None, background=False, logdir=None, log_prefix=''):
    import subprocess
    if not environment:
        environment = os.environ.copy()
    if logdir:
        logfile = os.path.join(logdir, log_prefix + '_$$.log')
        cmd += '>', logfile
    print("Working directory: {}".format(os.getcwd()))
    print("Command: {}".format(' '.join(cmd)))
    sys.stdout.flush()
    if background:
        process = subprocess.Popen(' '.join(cmd), env=environment, start_new_session=background, shell=True,
                                   stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        process = subprocess.Popen(cmd, env=environment,)

    print("Process ID: {}".format(process.pid))
    print()
    sys.stdout.flush()

    if not background:
        process.wait()


def pose_get_chain(pose, chain):
    from pyrosetta.rosetta.core.pose import append_pose_to_pose
    from pyrosetta.rosetta.core.pose import remove_nonprotein_residues
    '''
    Given a Rosetta pose and a chain letter, get only the relevant chain
    pose.
    '''
    remove_nonprotein_residues(pose)
    poses = []
    for i in range(1, pose.num_chains() + 1):
        chainpose = pose.split_by_chain(i)
        info = chainpose.pdb_info().pose2pdb(1)
        if info.split(' ')[1].strip() in chain:
            # if chainpose.size() < 5:
                # raise('Error: chain {} too small.'.format(chain))
            # else:
            poses.append(chainpose)
    pose = poses[0]
    if len(poses) > 1:
        for chainpose in poses[1:]:
            append_pose_to_pose(pose, chainpose)

    return pose


def trim_benchmark_df(df, col='name', col2='target'):
    # final_names = ['1d9c', '1e7l', '1i36', '1kqp', '1mty', '1nh2', '1p6x', '1um0', '1wui', '1xg0', '1yg6', '1ykd', '2dh4', '2e52', '2f1k', '2f4m', '2ftx', '2hp3', '2hzl', '2j91', '2ov9', '2pa8', '2q73', '2qib', '2qtq', '2ycl', '2yf3', '2yxh', '2zal', '3bpj', '3dhi', '3g5o', '3ilx', '3kkb', '3kra', '3sho', '3tos', '4g92', '4i0x', '4ics', '4lfg']
    final_names = ['1d9c', '1e7l', '1i36', '1kqp', '1mty', '1nh2', '1p6x', '1um0', '1yg6', '1ykd', '2dh4', '2e52', '2f4m', '2ftx', '2hp3', '2hzl', '2j91', '2ov9', '2pa8', '2q73', '2qib', '2qtq', '2ycl', '2yf3', '2yxh', '3bpj', '3dhi', '3g5o', '3ilx', '3kkb', '3kra', '4g92', '4i0x', '4lfg']
    df = df[df[col].isin(final_names)]
    df = df[~((df[col]=='1d9c') & (df[col2]=='B'))]
    return df


if __name__=='__main__':
    init()
    test_correct_resnums()
