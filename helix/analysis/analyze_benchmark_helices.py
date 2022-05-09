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
    start = min(row['rosetta_resis'])
    stop = max(row['rosetta_resis'])
    poses = []
    for i in range(start, stop - length + 1):
        this_stop = i + length
        this_target = target.clone()
        this_pose = rosetta.core.pose.append_subpose_to_pose(this_target, pose_clone, i, this_stop)
        poses.append(this_pose)

    return poses


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
            poses = make_bench_helix_pose(row['name'], row, length)
            for pose in poses:
                row = analyze_structures.analyze_pose(pose)
                row['length'] = f'len_{length}'
                print(row)
                outrows.append(row)

    outfile = os.path.join(output_folder, 'benchmark_scored.pkl')
    print(f'Saving to file: {outifle}')
    pd.DataFrame(outrows).to_pickle(outfile)