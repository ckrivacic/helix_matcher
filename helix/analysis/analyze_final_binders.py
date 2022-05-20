'''
Used to analyze the final binders from the 2022 interface design paper by Longxing Cao et al.

Usage:
    analyze_benchmark_interfaces.py <workspace> <output_folder> [options]

Options:
    --ntasks=INT  How many tasks total?  [default: 1]
    --task=INT  Run this task
'''
from pyrosetta import rosetta
# from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta import init
from pyrosetta import pose_from_file
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
from pyrosetta.rosetta.core.chemical import VariantType


def make_bench_helix_pose(pose, row):
    pose_clone = pose.clone()
    print(row)
    target = utils.pose_get_chain(pose, row['target'])
    print(target)
    # print('Target residue before', flush=True)
    # print(target.residue(1), flush=True)
    # rosetta.core.pose.remove_lower_terminus_type_from_pose_residue(target, 1)
    # print(target.residue(1), flush=True)
    chA = utils.pose_get_chain(pose, row['chain'])
    # print(chA.residue(chA.size()), flush=True)
    # rosetta.core.pose.remove_upper_terminus_type_from_pose_residue(chA, chA.size())
    # print(chA.residue(chA.size()), flush=True)
    print('TARGET POSE: ', target)
    print('chA POSE: ', chA)
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
    elif args['--task']:
        task = int(args['--task']) - 1
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
    inputs = sorted(glob.glob(
        os.path.join(output_folder, '*.pdb')))
    print('INPUTS')
    print(inputs)
    # benchmark_df = utils.safe_load(os.path.expanduser('~/software/helix_matcher/helix/benchmark/interface_finder/final_consolidated.pkl'))
    # groups = benchmark_df.groupby(['name','target','chain'])
    # names = sorted(list(groups.groups.keys()))
    # interval =  1 + (len(names) // int(args['--ntasks']))
    # start = task * interval
    # stop = task * interval + interval
    # names = names[start:stop]
    my_job = [inputs[task]]
    print('My job')
    print(my_job)

    outrows = []
    for pdb_file in my_job:
        initial_pose = pose_from_file(pdb_file)
        clone = initial_pose.clone()
        interface = score_pdb.PDBInterface(initial_pose, minimize=True, cst=True, is_pose=True)
        minimized_pose = interface.pose
        ref = create_score_function('ref2015')
        outrow = apply_filters(workspace, minimized_pose)
        outrow['total_score'] = ref(minimized_pose)
        outrow['chainA_size'] = minimized_pose.chain_end(1)
        outrow['name'] = f"{benchrow['name']}_{benchrow['target']}"
        outrow['sc_cst'] = False
        outrow['minimized']=True
        outrows.append(outrow)

        outrow_2 = apply_filters(workspace, clone)
        outrow_2['total_score'] = ref(minimized_pose)
        outrow_2['chainA_size'] = minimized_pose.chain_end(1)
        outrow_2['name'] = f"{benchrow['name']}_{benchrow['target']}"
        outrow_2['sc_cst'] = False
        outrow_2['minimized']=False
        outrows.append(outrow_2)

    outfile = os.path.join(output_folder, f'benchmark_scored_{task}.pkl')
    print(f'Saving to file: {outfile}')
    pd.DataFrame(outrows).to_pickle(outfile)


if __name__=='__main__':
    main()