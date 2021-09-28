from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta import pose_from_pdb
from pyrosetta import create_score_function
import pyrosetta
from distutils.dir_util import copy_tree
import docopt
import os, sys, glob
from helix import workspace as ws
from helix.utils import utils
from helix import big_jobs


def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector


def main():
    args = docopt.docopt(__doc__)
    workspace, job_info = big_jobs.initiate()
    if not hasattr(workspace, 'docking_directory'):
        raise Exception("Error: run_sge.py requires RIFWorkspaces as input. You "\
                "may have provided the root directory to this script "\
                "somehow.")

    total_jobs = len(glob.glob(workspace.focus_dir + '/patch_*/len_*'))
    print('TOTAL JOBS: {}'.format(total_jobs))
    num_tasks = total_jobs

    task_id = job_info['task_id']
    if not task_id:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0

    inputs = sorted(glob.glob(
        os.path.join(workspace.focus_dir, 'patch_*',
            workspace.scaffold_prefix + '*', 'docked_full', '*.pdb.gz')
        ))

    print('TASK: {}'.format(task_id))
    pdb = inputs[task_id]
    target = workspace.target_path_clean
    os.system('ls ???_????_*_*.pdb > input_list')
    # inputs = []
    # with open('input_list', 'r') as f:
        # for line in f:
            # inputs.append(line.strip())
    # for pdb in inputs:
    pose = pose_from_pdb(pdb)
    ref = create_score_function('ref2015')
    fastdes = pyrosetta.rosetta.protocols.denovo_designs.movers.FastDesign(ref)

    selector = residue_selector.ChainSelector('B')
    not_selector = residue_selector.NotResidueSelector(selector)
    tf = TaskFactory()
    no_packing = operation.PreventRepackingRLT()
    static = operation.OperateOnResidueSubset(no_packing,
            not_selector)
    notaa = operation.ProhibitSpecifiedBaseResidueTypes(
            strlist_to_vector1_str(['GLY']),
            selector)
    tf.push_back(static)
    tf.push_back(notaa)
    packertask = tf.create_task_and_apply_taskoperations(pose)
    print('REPACK')
    print(packertask.repacking_residues())
    print('DESIGN')
    print(packertask.designing_residues())

    fastdes.set_task_factory(tf)
    fastdes.set_scorefxn(ref)
    fastdes.apply(pose)
    score = ref(pose)
    pose.dump_pdb(pdb)

    exe = os.path.join(
            workspace.rosetta_dir, 'source', 'bin',
            'FlexPepDocking.{}'.format(suffix)
            )
    # if not os.path.exists('docked_full/'):
        # os.path.makedirs('docked_full', exist_ok=True)
    folder = os.path.dirname(
            os.path.abspath(pdb)
            )
    os.chdir(folder)

    cmd = [exe, '-in:file:s', pdb, '-scorefile',
            'score.sc',
            '-out:pdb_gz', '-lowres_preoptimize',
            '-flexPepDocking:pep_refine',
            '-flexPepDocking:flexpep_score_only', '-ex1', '-ex2aro',
            '-use_input_sc', 'unboundrot', target]
            # '-out:prefix', 'docked_full/',
    utils.run_command(cmd)
