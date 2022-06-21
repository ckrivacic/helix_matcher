"""
Script to analyze forward-folded designs (which are in silent file format).
Needs a workspace to get a few paths.

Usage:
    analyze_structures.py <workspace> <roseasy_focusdir> [options]

Options:
    --task=INT  Only run a certain task
    --designs-per-task=INT  How many designs to analyze per task  [default: 5]
    --fsf  Run the fragmetn score filter
"""

import numpy as np
from pyrosetta.rosetta.core.import_pose.pose_stream import SilentFilePoseInputStream
from pyrosetta.rosetta.core.scoring import calpha_superimpose_pose
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
from pyrosetta import init
from pyrosetta import Pose
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
import docopt
import os, glob
import pandas as pd
from helix import workspace as ws
from helix.design.design_interface import run_monomer_filters
from helix.design.design_interface import calculate_fsf
from roseasy import pipeline


def main():
    args = docopt.docopt(__doc__)

    if args['--task']:
        task_id = int(args['--task']) - 1
    elif 'SGE_TASK_ID' in os.environ:
        task_id = int(os.environ['SGE_TASK_ID']) - 1
    else:
        task_id = 0
    workspace = ws.workspace_from_dir(args['<workspace>'])
    roseasy_workspace = pipeline.workspace_from_dir(args['<roseasy_focusdir>'])
    dalphaball = os.path.join(workspace.rosetta_dir,
                              'source', 'external', 'DAlpahBall',
                              'DAlphaBall.gcc')
    ss_vall = workspace.find_path('ss_grouped_vall_all.h5')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
         ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false ' \
         '-indexed_structure_store:fragment_store {}'.format(dalphaball, ss_vall))

    inputs = sorted(glob.glob(roseasy_workspace.output_dir + '/*/*.out'))
    total_jobs = len(inputs)
    print('TOTAL JOBS: {}'.format(total_jobs), flush=True)
    tasks = int(args['--designs-per-task'])
    my_inputs = inputs[task_id*tasks:(task_id+1)*tasks]

    df = pd.DataFrame()
    for input in my_inputs:
        rowlist = []
        output_folder = os.path.dirname(input)
        pose = Pose()
        print('Loading input:')
        print(input)
        pis = SilentFilePoseInputStream(input)
        target = os.path.basename(os.path.dirname(input)).split('_')[0]
        input_name = os.path.basename(os.path.dirname(input))
        input_file = os.path.join(roseasy_workspace.input_dir, input_name + '.pdb.gz')
        init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
             ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false ' \
             '-indexed_structure_store:fragment_store {} '\
             '-s {}'.format(dalphaball, ss_vall, input_file))
        input_pose = pose_from_file(input_file)
        i = 0
        while pis.has_another_pose():
            pis.fill_pose(pose)
            # pdb = inputs[input_idx]
            # row = analyze_pose(pose, chA='A', chB='B')
            # row = apply_filters(workspace, pose)
            row = {}
            sfxn = create_score_function('ref2015')
            row['total_score'] = sfxn(pose)
            row = run_monomer_filters(workspace, pose, row)
            row['file'] = inputs[task_id]
            row['file_idx'] = i
            row['target'] = target
            row['design'] = input_file
            row['design_relpath'] = os.path.relpath(row['design'], start=roseasy_workspace.root_dir)
            row['design_name'] = '_'.join(input_name.split('_')[1:])
            model_no_idx = row['design_name'].split('_').index('model') + 1
            row['model_number'] = row['design_name'].split('_')[model_no_idx]
            row['chainA_size'] = pose.size()
            # row['descriptor'] = pis.get_last_pose_descriptor_string()
            # row['target'] = row['descriptor'].split('/')[0].split('_')[-1]
            row['silent_file'] = input
            row = run_more_filters(row, pose, input_pose, roseasy_workspace, input_name, task_id, args['--fsf'])
            rowlist.append(row)
            pis.next_struct()
            i += 1

        df = pd.concat([df, pd.DataFrame(rowlist)], ignore_index=True)
    pickle_outpath = os.path.join(output_folder, f'{input_name}_{task_id}.pkl')
    print(df, flush=True)
    print('Saving as {}'.format(pickle_outpath), flush=True)
    # if not os.path.exists(pickle_outdir):
    #     os.makedirs(pickle_outdir, exist_ok=True)
    # dataframe_out = os.path.join(pickle_outdir, f'{input_name}_{task_id}.pkl')
    df.to_pickle(pickle_outpath)


def get_json(pdb_path, workspace):
    import json
    model_no = os.path.basename(pdb_path).split('_')[2]
    lhl_folder = os.path.join(workspace.root_dir, '..', 'regenerated_data_sets_2020_03',
                              'sequence_design_for_LHL_reshaping_2lv8_two_LHL',
                              'selected_designs_for_state_count')
    insertion_file = os.path.join(lhl_folder, f'insertion_points_{model_no}.json')
    with open(insertion_file, 'r') as f:
        insertions = json.load(f)
    return insertions


def run_more_filters(row, pose, input_pose, workspace, filename, taskid, fsf=False):
    sup_rmsd = calpha_superimpose_pose(input_pose, pose)
    ca_rmsd = CA_rmsd(pose, input_pose)
    aa_rmsd = all_atom_rmsd(pose, input_pose)
    row['sup_rmsd'] = sup_rmsd
    row['ca_rmsd'] = ca_rmsd
    row['all_atom_rmsd'] = aa_rmsd
    i = 0
    test_run = False
    if fsf:
        for insertion in get_json(filename, workspace):
            i += 1
            try:
                row[f"frag_score_filter_{i}"] = calculate_fsf(workspace, pose, insertion,
                                                              f"{filename}_{str(taskid)}_{i}",
                                                              test_run=test_run)
                # test_run=True)
            except:
                row[f"frag_score_filter_{i}"] = np.nan

    return row


if __name__=='__main__':
    main()