"""
Script to analyze forward-folded designs (which are in silent file format).
Needs a workspace to get a few paths.

Usage:
    analyze_structures.py <workspace> <roseasy_focusdir> [options]

Options:
    --task=INT  Only run a certain task
    --designs-per-task=INT  How many designs to analyze per task  [default: 5]
"""

from pyrosetta.rosetta.core.import_pose.pose_stream import SilentFilePoseInputStream
from pyrosetta import init
from pyrosetta import Pose
import docopt
import os, glob
import pandas as pd
from helix import workspace as ws
from helix.design.design_interface import run_monomer_filters
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
        pose = Pose()
        print('Loading input:')
        print(input)
        pis = SilentFilePoseInputStream(input)
        target = os.path.basename(os.path.dirname(input)).split('_')[0]
        i = 0
        while pis.has_another_pose():
            pis.fill_pose(pose)
            # pdb = inputs[input_idx]
            if pose.num_chains() < 2:
                continue
            # row = analyze_pose(pose, chA='A', chB='B')
            # row = apply_filters(workspace, pose)
            row = {}
            row = run_monomer_filters(workspace, pose, row)
            row['file'] = inputs[task_id]
            row['file_idx'] = i
            row['target'] = args['<folder>']
            # row['descriptor'] = pis.get_last_pose_descriptor_string()
            # row['target'] = row['descriptor'].split('/')[0].split('_')[-1]
            row['silent_file'] = inputs[task_id]
            rowlist.append(row)
            pis.next_struct()
            i += 1

        df = pd.concat([df, pd.DataFrame(rowlist)], ignore_index=True)
    pickle_outdir = os.path.join(roseasy_workspace.focus_dir, 'analysis')
    print(df, flush=True)
    print('Saving in folder {}'.format(pickle_outdir), flush=True)
    if not os.path.exists(pickle_outdir):
        os.makedirs(pickle_outdir, exist_ok=True)
    dataframe_out = os.path.join(pickle_outdir, f'{target}_{task_id}.pkl')
    df.to_pickle(dataframe_out)


if __name__=='__main__':
    main()