'''
Make docking benchmark plots.

Usage:
    helix docking_benchmark <patchman_workspace> <rifdock_workspace> [options]


Options:

'''

from helix.utils import utils
import docopt
import os
import helix.workspace as ws


def plot_distribution(patch, rif, benchmark, args):
    '''Plot full distributions for the patchman and rifdock datasets'''
    for name, group in benchmark.groupby(['name', 'target']):
        patch_group = patch[(patch['name'] == name[0]) & (patch['target'] == name[1])]
        rif_group = rif[(rif['name'] == name[0]) & (rif['target'] == name[1])]
        for idx, row in group.iterrows():
            benchmark_resis = row['benchmark_resis']
            patch_subgroup = patch_group[patch_group['benchmark_resis'] == benchmark_resis]
            rif_subgroup = rif_group[rif_group['benchmark_resis'] == benchmark_resis]
            best_patchman = patch_subgroup.sort_values(by='best_rmsd', ascending=True).iloc[0]
            best_rifdock = rif_subgroup.sort_values(by='best_rmsd', ascending=True).iloc[0]
            print(best_patchman)
            print(best_rifdock)


def get_benchmark_resis(row):
    rosetta_resis = row['rosetta_resis']
    start = min(rosetta_resis)
    stop = max(rosetta_resis)
    return (start, stop)


def main():
    args = docopt.docopt(__doc__)
    bench_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
                 '..', 'benchmark', 'interface_finder', 'final_consolidated.pkl'
         )
    benchmark = utils.safe_load(bench_path)
    patchman_workspace = ws.workspace_from_dir(args['<patchman_workspace>'])
    rifdock_workspace = ws.workspace_from_dir(args['<rifdock_workspace>'])

    patchman_df = utils.safe_load(os.path.join(
        patchman_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final.pkl'
    ))
    rifdock_df = utils.safe_load(os.path.join(
        rifdock_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final.pkl'
    ))

    plot_distribution(patchman_df, rifdock_df, benchmark, args)