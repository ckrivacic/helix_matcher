'''
Make docking benchmark plots.

Usage:
    helix docking_benchmark <patchman_workspace> <rifdock_workspace> [options]


Options:
    --length=INT, -l  Length of helices to be evaluated. Can be either 14, 28, or 0. 0 means all lengths are
    evaluated together.  [default: 0]
    --trim=INT, -t  Trim the dataframe such that only benchmark targets that have at least two helices greater than
    the provided length are analyzed.

'''

from helix.utils import utils
import docopt
import os
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import helix.workspace as ws


def get_best_rmsds(patch, rif, benchmark, args):
    '''get the best RMSD for each benchmark helix'''
    no_patch_match = []
    no_rif_match = []
    outrows = []
    if args['--length'] == '14':
        patch = patch[patch['patch_len'] == 'len_14']
        rif = rif[rif['patch_len'] == 'len_4turn_dock_helix']
    elif args['--length'] == '28':
        patch = patch[patch['patch_len'] == 'len_28']
        rif = rif[rif['patch_len'] == 'len_8turn_dock_helix']
    for name, group in benchmark.groupby(['name', 'target']):
        patch_group = patch[(patch['name'] == name[0]) & (patch['target'] == name[1])]
        rif_group = rif[(rif['name'] == name[0]) & (rif['target'] == name[1])]
        for idx, row in group.iterrows():
            benchmark_resis = row['start_stop']
            patch_subgroup = patch_group[patch_group['start_stop'] == benchmark_resis]
            rif_subgroup = rif_group[rif_group['start_stop'] == benchmark_resis]
            best_patchman = patch_subgroup.sort_values(by='best_rmsd', ascending=True)
            best_rifdock = rif_subgroup.sort_values(by='best_rmsd', ascending=True)
            if best_patchman.shape[0] == 0 or best_rifdock.shape[0] == 0:
                if best_patchman.shape[0] == 0:
                    no_patch_match.append(row['name'])
                if best_rifdock.shape[0] == 0:
                    no_rif_match.append(row['name'])
            else:
                best_patchman = best_patchman.iloc[0]
                best_rifdock = best_rifdock.iloc[0]
                outrows.append({'name': row['name'], 'target': row.target, 'start_stop': row.start_stop,
                                      'rmsd': best_patchman.best_rmsd, 'protocol': 'PatchMAN'})
                outrows.append({'name': row['name'], 'target': row.target, 'start_stop': row.start_stop,
                                'rmsd': best_rifdock.best_rmsd, 'protocol': 'RIFDock'})

    print(no_patch_match)
    print(no_rif_match)
    outdf = pd.DataFrame(outrows)
    return outdf


def plot_distribution(df, args):
    '''Plot the distribution for pathcman and rifdock'''
    sns.histplot(data=df, x='rmsd', hue='protocol', stat='probability', common_norm=False)
    plt.show()


def get_benchmark_resis(row):
    rosetta_resis = row['rosetta_resis']
    start = min(rosetta_resis)
    stop = max(rosetta_resis)
    return (start, stop)


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)

    outpath = 'benchmark_data_{}.pkl'.format(args['--length'])
    bench_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..', 'benchmark', 'interface_finder', 'final_consolidated.pkl'
    )
    if not os.path.exists(outpath) or args['--trim']:
        benchmark = utils.safe_load(bench_path)
    if not os.path.exists(outpath):
        benchmark['start_stop'] = benchmark.apply(get_benchmark_resis, axis=1)
        patchman_workspace = ws.workspace_from_dir(args['<patchman_workspace>'])
        rifdock_workspace = ws.workspace_from_dir(args['<rifdock_workspace>'])

        patchman_df = utils.safe_load(os.path.join(
            patchman_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final.pkl'
        ))
        patchman_df['start_stop'] = patchman_df.apply(get_benchmark_resis, axis=1)
        rifdock_df = utils.safe_load(os.path.join(
            rifdock_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final.pkl'
        ))
        rifdock_df['start_stop'] = rifdock_df.apply(get_benchmark_resis, axis=1)
        df = get_best_rmsds(patchman_df, rifdock_df, benchmark, args)
        df.to_pickle(outpath)
        print(df.shape)
    else:
        df = utils.safe_load(outpath)

    if args['--trim']:
        benchmark['helix_length'] = benchmark.apply(lambda x: len(x['pdb_resis']), axis=1)

        def count(benchdf, length=14):
            return benchdf[benchdf['helix_length'] > length].shape[0]

        trim = []
        trim_cutoff = int(args['--trim'])
        for name, group in benchmark.groupby(['name', 'target']):
            if count(group, length=trim_cutoff) < 2:
                trim.append(name)

        print('The following targets will be trimmed:')
        print(trim)
        df['tup'] = df.apply(lambda x: (x['name'], x['target']), axis=1)
        for name in trim:
            print(df[df['tup'] == name])

        print(df.shape)
        # df = df[(~df['name'].isin(trim_name)) & (~df['target'].isin(trim_chain))]
        df = df[~df['tup'].isin(trim)]
        print(df.shape)

    plot_distribution(df, args)