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
        print(name)
        # patch_group = patch[(patch['name'] == )]


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

    plot_distribution(patch, rif, benchmark, args)