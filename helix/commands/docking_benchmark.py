'''
Make docking benchmark plots.

Usage:
    helix docking_benchmark <patchman_workspace> <rifdock_workspace> [options]


Options:
    --length=INT, -l  Length of helices to be evaluated. Can be either 14, 28, or 0. 0 means all lengths are
    evaluated together.  [default: 0]
    --trim=INT, -t  Trim the dataframe such that only benchmark targets that have at least two helices greater than
    the provided length are analyzed.
    --subangstrom, -s  Plot the percent sub-angstrom, also showing decreasing % ID
    --patch-id, -p  Plot percent sub-Angstrom based on patch % ID
'''

from helix.utils import utils
from helix.utils import homology
import docopt
import os
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import matplotlib.pyplot as plt
import helix.workspace as ws
from helix.utils.colors import palette


def get_patchman_pdbid(row):
    return os.path.basename(row['patchman_file']).split('_')[1].upper()


def calc_patch_percent_id(row):
    seqlen = 0
    identity = 0
    for idx, aa in enumerate(row.match_sequence):
        seqlen += 1
        if row.patch_sequence[idx] == aa:
            identity += 1
    if seqlen == 0:
        print('No designed sequences?')
        print('FILE: {}'.format(row.design_file))
        print('DESIGNED RESIS: {}'.format(row.designed_residues))
        print('HELIX RESIS: {}'.format(row.helix_resis))
        return np.nan

    return 100 * identity / seqlen


def get_best_rmsds(patch, rif, benchmark, args, cutoff=False, patch_id=False):
    '''get the best RMSD for each benchmark helix'''
    no_patch_match = {}
    no_rif_match = []
    outrows = []
    if args['--length'] == '14':
        patch = patch[patch['patch_len'] == 'len_14']
        rif = rif[rif['patch_len'] == 'len_4turn_dock_helix']
    elif args['--length'] == '28':
        patch = patch[patch['patch_len'] == 'len_28']
        rif = rif[rif['patch_len'] == 'len_8turn_dock_helix']
    for name, group in benchmark.groupby(['name', 'target']):
        print(name)
        patch_group = patch[(patch['name'] == name[0]) & (patch['target'] == name[1])]
        rif_group = rif[(rif['name'] == name[0]) & (rif['target'] == name[1])]
        print(rif_group)
        for idx, row in group.iterrows():
            benchmark_resis = row['start_stop']
            patch_subgroup = patch_group[patch_group['start_stop'] == benchmark_resis]
            if cutoff:
                cutoffs = set(patch_subgroup['cutoff'])
                column = 'cutoff'
            elif patch_id:
                cutoffs = [10, 30, 50, 70, 90, 95, 100, 110]
                column = 'patch_percent_id'
            else:
                cutoffs = [0]
            for c in cutoffs:
                no_patch_match[c] = []
                if cutoff:
                    patch_subgroup = patch_group[patch_group[column] == c]
                if patch_id:
                    patch_subgroup = patch_group[patch_group[column] < c]
                best_patchman = patch_subgroup.sort_values(by='best_rmsd', ascending=True)
                if best_patchman.shape[0] == 0:
                    no_patch_match[c].append(row['name'])
                else:
                    best_patchman = best_patchman.iloc[0]
                    outrows.append({'name': row['name'], 'target': row.target, 'start_stop': row.start_stop,
                                    'rmsd': best_patchman.best_rmsd, 'protocol': 'PatchMAN', 'cutoff': c})

            rif_subgroup = rif_group[rif_group['start_stop'] == benchmark_resis]
            best_rifdock = rif_subgroup.sort_values(by='best_rmsd', ascending=True)
            if best_rifdock.shape[0] == 0:
                no_rif_match.append(row['name'])
            else:
                best_rifdock = best_rifdock.iloc[0]
                outrows.append({'name': row['name'], 'target': row.target, 'start_stop': row.start_stop,
                                'rmsd': best_rifdock.best_rmsd, 'protocol': 'RIFDock', 'cutoff': 0})

    print('No patch match:')
    print(no_patch_match)
    print('No RIF match:')
    print(no_rif_match)
    outdf = pd.DataFrame(outrows)
    return outdf


def plot_distribution(df, args):
    '''Plot the distribution for pathcman and rifdock'''
    sns.histplot(data=df, x='rmsd', hue='protocol', stat='probability', common_norm=False)
    plt.show()


def plot_subA(df, args):
    '''Plot percentage of sub-Angstrom benchmark helices'''
    df = df.drop_duplicates(ignore_index=True)
    order = ['RIFDock_0', 'PatchMAN_110', 'PatchMAN_100', 'PatchMAN_95', 'PatchMAN_90', 'PatchMAN_70', 'PatchMAN_50',
             'PatchMAN_30', 'PatchMAN_10']
    labels = ['RIFDock', 'PatchMAN', '100%', '95%', '90%', '70%', '50%', '30%', '10%']
    colors = [palette['red'], palette['blue'], '#4395c9', '#58a1cf', '#6dadd5', '#82b8db', '#97c4e1', '#abd0e7', '#c0dbed']
    sns.set_palette(colors)

    subA_df = []
    for name, group in df.groupby(['protocol', 'cutoff']):
        protocol = f'{name[0]}_{name[1]}'
        fraction_subA = group[group['rmsd'] < 1.0].shape[0] / group.shape[0]
        subA_df.append({
            'protocol': protocol,
            'fraction_subA': fraction_subA
        })
    # plt.figure(figsize=(3,3), dpi=300)
    fig, ax = plt.subplots(figsize=(3,3), dpi=300)
    subA_df = pd.DataFrame(subA_df)
    sns.barplot(x='protocol', y='fraction_subA', data=subA_df, order=order)
    ax.set_xticklabels(labels)
    plt.tight_layout()
    plt.xticks(rotation=70)
    plt.xlabel(None)
    plt.ylabel("Fraction sub-Å")
    plt.show()


def get_benchmark_resis(row):
    rosetta_resis = row['rosetta_resis']
    start = min(rosetta_resis)
    stop = max(rosetta_resis)
    return (start, stop)


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)

    outpath = 'benchmark_data_{}'.format(args['--length'])
    if args['--subangstrom']:
        outpath += '_homology'
    if args['--patch-id']:
        outpath += '_patch_ID'
    if args['--trim']:
        outpath += '_trimmed'
    outpath += '.pkl'
    bench_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..', 'benchmark', 'interface_finder', 'final_consolidated.pkl'
    )
    patchman_workspace = ws.workspace_from_dir(args['<patchman_workspace>'])
    rifdock_workspace = ws.workspace_from_dir(args['<rifdock_workspace>'])
    if not os.path.exists(outpath) or args['--trim']:
        benchmark = utils.safe_load(bench_path)
    if not os.path.exists(outpath):
        benchmark['start_stop'] = benchmark.apply(get_benchmark_resis, axis=1)

        print('Loading patchman DF')
        patchman_df = utils.safe_load(os.path.join(
            patchman_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final_aligned.pkl'
        ))
        patchman_df['start_stop'] = patchman_df.apply(get_benchmark_resis, axis=1)
        print('Loading RIFDock DF')
        rifdock_df = utils.safe_load(os.path.join(
            rifdock_workspace.root_dir, 'rifdock_outputs', 'benchmark_results_reverse', 'final.pkl'
        ))
        rifdock_df['start_stop'] = rifdock_df.apply(get_benchmark_resis, axis=1)
        if args['--patch-id']:
            print('Calculating patch/match percent ID...')
            patchman_df['patch_percent_id'] = patchman_df.apply(calc_patch_percent_id, axis=1)
        if not args['--subangstrom']:
            print('Finding best RMSD for each benchmark helix')
            df = get_best_rmsds(patchman_df, rifdock_df, benchmark, args, patch_id=args['--patch-id'])
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

    if args['--subangstrom']:
        if not os.path.exists(outpath):
            homology_outpath = os.path.join(
                patchman_workspace.rifdock_outdir,
                'benchmark_results_reverse',
                'patchman_results_homology.pkl'
            )
            if not os.path.exists(homology_outpath):
                cutoffs = [10, 30, 50, 70, 80, 90, 95, 100, 110]
                print('Figuring out PDBIDs...')
                patchman_df['patchman_pdbid'] = patchman_df.apply(get_patchman_pdbid, axis=1)
                homology_dfs = []
                for cutoff in cutoffs:
                    print('Finding homologs within {}% cutoff'.format(cutoff))
                    for name, group in patchman_df.groupby(['name', 'target']):
                        print('Finding homologs for {} with {} percent ID'.format(name, cutoff))
                        pdbid = name[0]
                        chain = name[1]
                        homologs = homology.find_homologs(pdbid, cutoff,
                                                          chain=chain)
                        group = group[~group['patchman_pdbid'].isin(homologs)]
                        group['cutoff'] = cutoff
                        homology_dfs.append(group)
                homology_df = pd.concat(homology_dfs)
                homology_df.to_pickle(homology_outpath)
            else:
                print('Loading homology df...')
                homology_df = utils.safe_load(homology_outpath)
            print(homology_df)
            rmsd_df = get_best_rmsds(homology_df, rifdock_df, benchmark, args, cutoff=args['--subangstrom'], patch_id=args['--patch-id'])
            rmsd_df.to_pickle(outpath)
        else:
            rmsd_df = utils.safe_load(outpath)
        plot_subA(rmsd_df, args)
    elif args['--patch-id']:
        plot_subA(df, args)

    else:
        plot_distribution(df, args)