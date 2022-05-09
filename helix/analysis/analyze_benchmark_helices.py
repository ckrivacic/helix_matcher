'''
Usage:
    analyze_benchmark_helices.py <output_folder> [options]
'''
from pyrosetta import rosetta
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
from pyrosetta.rosetta.core.scoring.dssp import Dssp
import docopt
import os, sys, glob
import pandas as pd
from helix.utils import utils
from helix.analysis import analyze_structures
from helix.rifdock import interface
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects


def make_bench_helix_pose(pose, row, length):
    pose_clone = pose.clone()
    target = utils.pose_get_chain(pose, row['target'])
    print('TARGET POSE: ', target)
    start = min(row['rosetta_resis'])
    stop = max(row['rosetta_resis'])
    if stop - start + 1 < length:
        length = stop - start
    poses = []
    print('RESIDUE RANGE: ', start, stop)
    for i in range(start, stop - length + 1):
        this_stop = i + length
        this_target = target.clone()
        print('THIS TARGET: ', this_target)
        print('TARGET SIZE: ', this_target.size())
        rosetta.core.pose.append_subpose_to_pose(this_target, pose_clone, i, this_stop, True)
        print('THIS POSE: ', this_target)
        print('THIS POSE SIZE: ', this_target.size())
        print('CHAIN 1 letter: ', this_target.pdb_info().pose2pdb(1))
        print('CHAIN 2 letter: ', this_target.pdb_info().pose2pdb(this_target.chain_begin(2)))
        print('CHAIN 2 SEQ: ', this_target.chain_sequence(2))
        poses.append(this_target)

    return poses, length


def main():
    args = docopt.docopt(__doc__)
    rosetta_dir = os.path.expanduser('~/software/rosetta/')
    dalphaball = os.path.join(rosetta_dir,
                              'source', 'external', 'DAlpahBall',
                              'DAlphaBall.gcc')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
         ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false'.format(dalphaball))
    output_folder = args['<output_folder>']
    benchmark_df = utils.safe_load(os.path.expanduser('~/software/helix_matcher/helix/benchmark/interface_finder/final_consolidated.pkl'))

    outrows = []
    for idx, row in benchmark_df.iterrows():
        print(row)
        for length in [14, 28]:
            initial_pose = utils.pose_from_wynton(row['name'])
            poses, length = make_bench_helix_pose(initial_pose, row, length)
            print('POSES', poses)
            for pose in poses:
                row = analyze_structures.analyze_pose(pose, row['target'], row['chain'])
                row['length'] = length
                print(row)
                outrows.append(row)

    outfile = os.path.join(output_folder, 'benchmark_scored.pkl')
    print(f'Saving to file: {outfile}')
    pd.DataFrame(outrows).to_pickle(outfile)


if __name__=='__main__':
    main()