'''
Usage:
    analyze_benchmark_interfaces.py <workspace> <output_folder> [options]

Options:
    --ntasks=INT  How many tasks total?  [default: 1]
'''
from pyrosetta import rosetta
# from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta import init
# from pyrosetta import pose_from_file
# from pyrosetta import create_score_function
# from pyrosetta.rosetta.core.scoring.dssp import Dssp
import docopt
import os, sys, glob
import pandas as pd
from helix.utils import utils
from helix.benchmark import score_pdb
# from helix.rifdock import interface
# from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta import create_score_function
import helix.workspace as ws
from helix.design.design_interface import apply_filters


def make_bench_helix_pose(pose, row):
    pose_clone = pose.clone()
    print(row)
    target = utils.pose_get_chain(pose, row['target'])
    chA = utils.pose_get_chain(pose, row['chain'])
    print('TARGET POSE: ', target)
    rosetta.core.pose.append_pose_to_pose(chA, target, True)
    for resnum in range(chA.chain_begin(1), chA.size() + 1):
        chA.pdb_info().chain(resnum, 'B')
    for resnum in range(1, chA.chain_end(1) + 1):
        chA.pdb_info().chain(resnum, 'A')

    return chA


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    if 'SGE_TASK_ID' in os.environ:
        task = int(os.environ['SGE_TASK_ID']) - 1
    else:
        task = 0

    # rosetta_dir = os.path.expanduser('~/rosetta/')
    # dalphaball = os.path.join(rosetta_dir,
    #                           'source', 'external', 'DAlpahBall',
    #                           'DAlphaBall.gcc')
    # init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
    #      ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false'.format(dalphaball))
    dalphaball = os.path.join(workspace.rosetta_dir,
                              'source', 'external', 'DAlpahBall',
                              'DAlphaBall.gcc')
    ss_vall = workspace.find_path('ss_grouped_vall_all.h5')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
         ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false ' \
         '-indexed_structure_store:fragment_store {}'.format(dalphaball, ss_vall))
    output_folder = args['<output_folder>']
    benchmark_df = utils.safe_load(os.path.expanduser('~/software/helix_matcher/helix/benchmark/interface_finder/final_consolidated.pkl'))
    interval =  1 + (benchmark_df.shape[0] // int(args['--ntasks']))
    start = task * interval
    stop = task * interval + interval
    benchmark_df = benchmark_df.iloc[start:stop]

    outrows = []
    for idx, benchrow in benchmark_df.iterrows():
        print(benchrow)
        initial_pose = utils.pose_from_wynton(benchrow['name'])
        pose = make_bench_helix_pose(initial_pose, benchrow)
        interface = score_pdb.PDBInterface(pose, minimize=True, cst=True, is_pose=True)
        minimized_pose = interface.pose
        ref = create_score_function('ref2015')
        outrow = apply_filters(workspace, minimized_pose)
        outrow['total_score'] = ref(minimized_pose)
        outrow['chainA_size'] = minimized_pose.chain_end(1)
        outrow['name'] = f"{benchrow['name']}_{benchrow['target']}"
        outrow['sc_cst'] = False
        outrows.append(outrow)

    outfile = os.path.join(output_folder, f'benchmark_scored_{task}.pkl')
    print(f'Saving to file: {outfile}')
    pd.DataFrame(outrows).to_pickle(outfile)


if __name__=='__main__':
    main()