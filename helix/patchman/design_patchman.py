"""
Script to run FastDesign + FlexPepDock on outputs of PatchMAN.

Usage:
    design_patchman.py <workspace> [options]

Options:
    --task=INT  Only run a specific task
    --delete  Delete target structures

"""

from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
import pyrosetta
from distutils.dir_util import copy_tree
import docopt
import os, sys, glob
import pandas as pd
import gzip
from helix import workspace as ws
from helix.utils import utils
from helix import big_jobs
from helix.rifdock import interface
# from helix.matching.scan_helices import contiguous_secstruct
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
import pymol


def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector


def main():
    dalphaball = os.path.join('/wynton', 'home', 'kortemme', 'krivacic',
            'rosetta', 'source', 'external', 'DAlpahBall',
            'DAlphaBall.gcc')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro'\
            ' -holes:dalphaball {}'.format(dalphaball))
    args = docopt.docopt(__doc__)
    try:
        workspace, job_info = big_jobs.initiate()
        # job_info['task_id'] = int(args['--task'])
        # job_info['inputs'] = sorted(glob.glob(
            # os.path.join(workspace.focus_dir, 'patch_*',
                # workspace.scaffold_prefix + '*', 'docked_full',
                # '*.pdb.gz'),
            # ))
    except:
        print('Maybe this is local?')
        workspace = ws.workspace_from_dir(args['<workspace>'])
        job_info = {
                'task_id': args['--task'],
                'inputs' : sorted(glob.glob(
                    os.path.join(workspace.focus_dir, 'patch_*',
                        workspace.scaffold_prefix + '*', 'docked_full',
                        '*.pdb.gz'),
            ))
                }

    if not hasattr(workspace, 'docking_directory'):
        raise Exception("Error: design_patchman.py requires RIFWorkspaces as input. You "\
                "may have provided the root directory to this script "\
                "somehow.")

    inputs = job_info['inputs']
    nstruct = 50
    total_jobs = len(inputs)
    print('TOTAL JOBS: {}'.format(total_jobs))

    # print('Job info')
    # print(job_info)
    if 'task_id' not in job_info:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0
    else:
        task_id = int(job_info['task_id'])
    start = task_id * nstruct
    stop = task_id * nstruct + nstruct

    target = pymol.cmd.load(workspace.target_path_clean)

    print('TASK: {}'.format(task_id))
    rowlist = []
    for input_idx in range(start, stop):
        pdb_save = os.path.relpath(inputs[input_idx],
                start=workspace.root_dir)
        pdb = os.path.abspath(inputs[input_idx])
        designed = False
        with gzip.open(pdb, 'rt') as f:
            for line in f:
                if line.startswith('pose'):
                    designed = True

        folder = os.path.dirname(
                os.path.abspath(pdb)
                )
        os.chdir(folder)
        target = workspace.target_path_clean
        # os.system('ls ???_????_*_*.pdb > input_list')
        # inputs = []
        # with open('input_list', 'r') as f:
            # for line in f:
                # inputs.append(line.strip())
        # for pdb in inputs:
        pose = pose_from_file(pdb)
        ref = create_score_function('ref2015')
        selector = residue_selector.ChainSelector('B')
        if not designed:
            fastdes = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign(ref)

            movemap = MoveMap()
            movemap.set_bb_true_range(pose.chain_begin(2),
                    pose.chain_end(2))
            movemap.set_chi_true_range(pose.chain_begin(2),
                    pose.chain_end(2))
            fastdes.set_movemap(movemap)

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

        # import platform
        # ostype = platform.system()
        # if ostype == 'Linux':
            # suffix = 'linuxgccrelease'
        # elif ostype == 'Darwin':
            # suffix = 'macosclangrelease'

        # exe = os.path.join(
                # workspace.rosetta_dir, 'source', 'bin',
                # 'FlexPepDocking.{}'.format(suffix)
                # )
        # if not os.path.exists('docked_full/'):
            # os.path.makedirs('docked_full', exist_ok=True)

        # cmd = [exe, '-in:file:s', pdb, '-scorefile',
                # 'score.sc',
                # '-out:pdb_gz', '-lowres_preoptimize',
                # '-flexPepDocking:pep_refine',
                # '-flexPepDocking:flexpep_score_only', '-ex1', '-ex2aro',
                # '-use_input_sc', '-unboundrot', target]
                # # '-out:prefix', 'docked_full/',
        # utils.run_command(cmd)

        pdb_basename = pdb.split('.')[0]
        # flexpep_file = pdb_basename + '_0001.pdb.gz'
        # flexpep_pose = pose_from_file(flexpep_file)
        flexpep_file = pdb
        flexpep_pose = pose
        chainB = pose.split_by_chain(2)
        pymol_mobile = pymol.cmd.load(pdb, 'mobile')
        pymol.cmd.align('mobile and chain A', 'target')
        pymol.cmd.save(pdb, 'mobile')
        pymol.cmd.delete('mobile')

        # Determine helical propensity
        ss_str = Dssp(chainB).get_dssp_secstruct()
        # secstruct = contiguous_secstruct(ss_str)
        percent_helical = ss_str.count('H') / len(ss_str)

        # Define filters
        buns_all = '''
        <BuriedUnsatHbonds 
            name="Buried Unsat [[-]]"
            report_all_heavy_atom_unsats="true" scorefxn="ref2015"
            cutoff="4" residue_surface_cutoff="20.0"
            ignore_surface_res="true" print_out_info_to_pdb="true"
            dalphaball_sasa="1" probe_radius="1.1" confidence="0"
            only_interface="true" />

        '''
        buns_sc = '''
        <BuriedUnsatHbonds 
            name="Buried Unsat Sidechains [[-]]"
            report_sc_heavy_atom_unsats="true" scorefxn="ref2015" cutoff="4"
            residue_surface_cutoff="20.0" ignore_surface_res="true"
            print_out_info_to_pdb="true"  dalphaball_sasa="1"
            probe_radius="1.1" confidence="0" only_interface="true" />

        '''
        npsa = '''
        <BuriedSurfaceArea name="Buried Nonpolar Surface Area [[+]]"
        select_only_FAMILYVW="true" filter_out_low="false"
        atom_mode="all_atoms"
        confidence="1.0"
        />
        '''
        packstat = '''
          <PackStat
            name="PackStat Score [[+]]"
            threshold="0"
          />
        '''
        buns_all_obj = XmlObjects.static_get_filter(buns_all)
        buns_sc_obj = XmlObjects.static_get_filter(buns_sc)
        npsa_obj = XmlObjects.static_get_filter(npsa)
        npsa_obj.set_residue_selector(selector)
        packstat_obj = XmlObjects.static_get_filter(packstat)
        buns_all_score = buns_all_obj.report_sm(flexpep_pose)
        buns_sc_score = buns_sc_obj.report_sm(flexpep_pose)
        npsa_score = npsa_obj.report_sm(flexpep_pose)
        packstat_score = packstat_obj.report_sm(flexpep_pose)

        score = ref(flexpep_pose)
        interface_scorer = interface.InterfaceScore(flexpep_pose)
        interface_score = interface_scorer.apply()
        n_hbonds = interface_scorer.n_hbonds
        row = {'patchman_file': pdb_save,
                'name': os.path.basename(flexpep_file),
                'size': flexpep_pose.size(),
                'pose_score': score,
                'interface_score': interface_score,
                'n_hbonds': n_hbonds,
                'buns_all': buns_all_score,
                'buns_sc': buns_sc_score,
                'npsa': npsa_score,
                'packstat': packstat_score,
                'percent_helical': percent_helical,
                }
        rowlist.append(row)

    df = pd.DataFrame(rowlist)
    print(df)
    pickle_outdir = os.path.join(workspace.focus_dir, 'scores')
    print('Saving in folder {}'.format(pickle_outdir))
    if not os.path.exists(pickle_outdir):
        os.makedirs(pickle_outdir, exist_ok=True)
    df.to_pickle(os.path.join(pickle_outdir,
        'task_{task}.pkl'.format(task=task_id)))

    # if args['--delete']:
        # os.remove(pdb)

if __name__=='__main__':
    main()
