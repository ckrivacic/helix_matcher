'''
Analyze design outputs
plot_type can be bar or scatter.

Usage:
    helix plot_designs <workspace> <plot_type> [options]

Options:
    --target=NAME  Only plot for a specific target
    --suffix=STR  Only plot for a specific suffix
    --compare-suffix  Plot different suffixes side by side or in different colors if scatterplot
    --xaxis=COLUMN, -x  Plot this column on the x-axis  [default: ca_rmsd]
    --yaxis=COLUMN, -y  Plot this column on the y-axis  [default: interface_dG]
    --hue=COLUMN  Color by values from this column
    --size=COLUMN  Resize scatterplot points by values from this column
    --additional-df  Compare to another dataframe (for ex., benchmark data or data from another design method)
    --force-reread, -f  Even if final.pkl exists in the design directory, reread all individual score files.
    --filter=FILEPATH  Filter the dataframe
    --keep-original  Keep original dataframe when filtering, plot both
    --compare-final  Compare to the "final" dataset from the Baker lab de novo project
    --compare-bench  Compare to analyzed benchmark itnerfaces
    --compare-design  Compare to design models from Baker interface project
    --perres-x  Calculate per residue metrics
    --perres-y  Calculate per residue metrics
    --model=NUMBER  Only plot this model #
'''
from helix.utils import utils
import helix.workspace as ws
from helix.utils import colors
import pandas as pd
from helix.utils import plotting
import docopt
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import os


def normalize(array, min_val=3, max_val=50):
    array = array - min_val
    x = max(array) / max_val
    array = array / x
    return array


global name_dict
name_dict = {
    '1DJS': 'FGFR2',
    '1MOX': 'EGFR',
    '1XIW': 'CD3δ',
    '2GY7': 'TIE2',
    '2IFG': 'TrkA',
    '3MJG_chX': 'PDGFR',
    '4O3V': 'VirB8',
    '4OGA': 'InsulinR',
    '5U8R': 'IGF1R',
    'N/A': '',
}
global folder_dict
folder_dict = {
    '1DJS': ['FGFR2_mb.pdb', 'FGFR2'],
    '1MOX': ['EGFRn_mb.pdb', 'EGFRc_mb.pdb', 'EGFR'],
    '1XIW': ['CD3d_mb.pdb', 'CD3', 'CD3d'],
    '2GY7': ['Tie2_mb.pdb', 'Tie2'],
    '2IFG': ['TrkA_mb.pdb', 'TrkA'],
    '3MJG_chX': ['PDGFR_mb.pdb', 'PDGFR'],
    '4O3V': ['VirB8_mb.pdb', 'VirB8'],
    '4OGA': ['InsulinR_mb.pdb', 'InsulinR', 'IR'],
    '5U8R': ['IGF1R_mb.pdb', 'IGF1R'],
    'N/A': '',
}


def parse_metric(metric_name):
    metric_dict = {
        'interface_dG': 'Interface ΔG (REU)',
        'contact_molecular_surface': 'CMS score',
        'n_hbonds': 'Hydrogen bonds'
    }
    if metric_name in metric_dict:
        return metric_dict[metric_name]
    else:
        name_list = '_'.split(metric_name)
        name_list[0] = name_list[0].capitalize()
        return ' '.join(name_list)


def violin_plot(dfs, workspace, args):
    '''Plot a violin plot; intention is to show distribution of interface Hbonds for designs with no BUNS'''
    from numpy import mean

    sns.categorical._Old_Violin = sns.categorical._ViolinPlotter
    class _My_ViolinPlotter(sns.categorical._Old_Violin):

        def __init__(self, *args, **kwargs):
            super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
            self.gray = 'black'
    sns.categorical._ViolinPlotter = _My_ViolinPlotter

    order = ['FGFR2', 'EGFR', 'CD3δ', 'TIE2', 'TrkA', 'PDGFR', 'VirB8', 'InsulinR', 'IGF1R']
    hue = 'Design source'
    fig, ax = plt.subplots(figsize=(4,3.5), dpi=300)
    def rename_target(row):
        if not args['--target'] == '1MOX':
            rename_dict = {
                'CD3': 'CD3δ',
                'CD3d': 'CD3δ',
                'Tie2': 'TIE2',
                'EGFRc': 'EGFR',
                'EGFRn': 'EGFR',
            }
            if row['target'] in rename_dict:
                return rename_dict[row['target']]
            else:
                return row['target']
        else:
            rename_dict = {''
                           'EGFRc': 'EGFRc',
                           'EGFRn': 'EGFRn',
                           }
            if row['target'] in rename_dict:
                return rename_dict[row['target']]
            else:
                return row['target']

    df = pd.concat(dfs, ignore_index=True)
    df['target'] = df.apply(lambda x: name_dict[x['design_file'].split('/')[1]], axis=1)
    df[hue] = 'HELIX'
    if args['--target'] == '1MOX':
        from copy import deepcopy
        print('Editing base dataframe')
        df = df[df['target'] == 'EGFR']
        df['target'] = 'EGFRc'
        df2 = deepcopy(df)
        df2['target'] = 'EGFRn'
        df = pd.concat([df, df2], ignore_index=True)
        print(df)
    if args['--compare-design']:
        design_df = utils.safe_load(os.path.expanduser(
            '~/intelligent_design/helix_workspaces/baker_design_models/final.pkl'
        ))
        design_df['target'] = design_df.apply(rename_target, axis=1)
        design_df[hue] = 'Miniproteins'
        if args['--filter']:
            from helix.utils.utils import parse_filter
            import yaml
            with open(args['--filter'], 'r') as f:
                filter = yaml.load(f.read())
                design_df = parse_filter(filter, design_df)
        if args['--target'] == '1MOX':
            design_df = design_df[(design_df['target'] == 'EGFRn') | (design_df['target'] == 'EGFRc')]
            print(design_df)
        df = pd.concat([df, design_df], ignore_index=True)

    if args['--compare-bench']:
        bench_df = utils.safe_load(os.path.expanduser(
            '~/intelligent_design/helix_workspaces/benchmark_interface_analysis/final.pkl'
        ))
        print('benchmark:')
        print(bench_df.columns)
        bench_df['true_name'] = bench_df.apply(lambda x: x['name'].split('_')[0], axis=1)
        bench_df = utils.trim_benchmark_df(bench_df, col='true_name')
        bench_df = bench_df[~bench_df['sc_cst']]
        bench_df[hue] = 'Natural interfaces'
        bench_df['target'] = 'IGF1R'
        df = pd.concat([df, bench_df])

    if args['--compare-final']:
        final_df = utils.safe_load(os.path.expanduser(
            '~/intelligent_design/helix_workspaces/design_models_final_combo_optimized/final.pkl'
        ))
        final_df = final_df[final_df['minimized']]
        print(final_df['name'].unique())
        final_df['target'] = final_df.apply(lambda x: x['name'].split('_')[0], axis=1)
        final_df['target'] = final_df.apply(rename_target, axis=1)
        final_df = final_df[final_df['target'].isin(df['target'].unique())]
    print(df)

    if args['--perres-x']:
        df[args['--xaxis']] = df.apply(
            lambda x: x[args['--xaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0,
            axis=1)
    if args['--perres-y']:
        df[args['--yaxis']] = df.apply(lambda x: x[args['--yaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)

    bw = 0.3
    palette = [colors.palette['teal'], colors.palette['yellow'], colors.palette['red']]
    sns_ax = sns.violinplot(data=df, x='target', y=args['--yaxis'],
                            order=order, hue=hue, bw=bw, dodge=True, split=True,
                            palette=palette, linewidth=1, inner=None)
    sns.pointplot(x='target', y=args['--yaxis'], data=df, estimator=mean, hue=hue, dodge=0.25, palette=['black'],
                  join=False, marker='hline', ci=0, scale=0, capsize=0.25, errwidth=1)
    sns.swarmplot(x='target', y=args['--yaxis'], data=final_df, palette=[colors.palette['purple']], marker='d',
                  edgecolor='black', linewidth=0.5)
                  #scale=0.5)

    sns_ax.set_xticklabels(order, ha='right', fontsize=8)
    plt.ylabel(parse_metric(args['--yaxis']), fontsize=8)
    plt.xticks(rotation=70)
    plt.legend([], [], frameon=False)
    ax.tick_params(axis='y', labelsize=6, length=2)
    ax.tick_params(axis='both', which='major', pad=2)
    plt.tight_layout()
    plt.show()


def scatterplot(dfs, workspace, args, use_matplotlib=True):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
    # '1b33_K_nativelike', '1b33_K_specialrot']
    if not args['--xaxis'] == 'protocol':
        order = None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    # print(df[args['--yaxis']].sort_values())

    if len(dfs) < 3:
        nrows = 1
        ncols = len(dfs)
    else:
        nrows = 3
        ncols = 3
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5, 3.5), dpi=300)
    axs = np.array(axs)
    # if args['--size']:
    #     size = args['--size']
    #     ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
    #                          hue=hue, size=size, sizes=(50,300), picker=True)
    for idx, ax in enumerate(axs.reshape(-1)):
        if idx > len(dfs)-1:
            continue
        df = dfs[idx]
        print('DF COLS')
        print(df.columns)
        if df.shape[0] > 0:
            target_name = df.iloc[0]['design_file'].split('/')[1]
        else:
            target_name = 'N/A'
        if args['--perres-x']:
            # df[args['--xaxis']] = df.apply(lambda x: x[args['--xaxis']] / x['chainA_size'], axis=1)
            df[args['--xaxis']] = df.apply(lambda x: x[args['--xaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)
        if args['--perres-y']:
            # df[args['--yaxis']] = df.apply(lambda x: x[args['--yaxis']] / x['chainA_size'], axis=1)
            df[args['--yaxis']] = df.apply(lambda x: x[args['--yaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)

        if args['--compare-design']:
            design_df = utils.safe_load(os.path.expanduser(
                '~/intelligent_design/helix_workspaces/baker_design_models/final.pkl'
            ))
            design_df = design_df[design_df['target'].isin(folder_dict[target_name])]
            if args['--filter']:
                from helix.utils.utils import parse_filter
                import yaml
                with open(args['--filter'], 'r') as f:
                    filter = yaml.load(f.read())
                    design_df = parse_filter(filter, design_df)
            if not design_df.empty:
                if args['--perres-x']:
                    design_df[args['--xaxis']] = design_df.apply(
                        lambda x: x[args['--xaxis']] / len(x['interfac_residues']), axis=1)
                if args['--perres-y']:
                    design_df[args['--yaxis']] = design_df.apply(
                        lambda x: x[args['--yaxis']] / len(x['interfac_residues']), axis=1)
                x = design_df[args['--xaxis']]
                y = design_df[args['--yaxis']]
                mpl_args = [x, y]
                kwargs = {}
                if args['--size']:
                    sizes = df[args['--size']]
                    sizes = normalize(sizes)
                    kwargs['s'] = sizes
                else:
                    kwargs['s'] = 1
                # if hue:
                kwargs['edgecolor'] = 'white'
                kwargs['linewidth'] = 0.05
                kwargs['c'] = colors.palette['yellow']
                points = ax.scatter(*mpl_args, picker=True, marker='s', **kwargs)

        if use_matplotlib:
            # Seaborn seems to reorder points, making the on_pick stuff useless. Back to MPL.
            x = df[args['--xaxis']]
            y = df[args['--yaxis']]
            mpl_args = [x, y]
            kwargs = {}
            if args['--size']:
                sizes = df[args['--size']]
                sizes = normalize(sizes)
                kwargs['s'] = sizes
            else:
                kwargs['s'] = 1
            if hue:
                kwargs['c'] = df[hue]
            else:
                kwargs['c'] = colors.palette['teal']
            kwargs['alpha'] = 0.5
            kwargs['edgecolor'] = 'white'
            kwargs['linewidth'] = 0.05
            points = ax.scatter(*mpl_args, picker=True, **kwargs)
            if args['--hue']:
                plt.colorbar(points, ax=ax)
            click = plotting.ClickablePlot(points, df, args, workspace)
        else:
            ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                                 hue=hue, picker=True, alpha=1.0)

        if args['--compare-bench']:
            bench_df = utils.safe_load(os.path.expanduser(
                '~/intelligent_design/helix_workspaces/benchmark_interface_analysis/final.pkl'
            ))
            print('benchmark:')
            print(bench_df.columns)
            bench_df['true_name'] = bench_df.apply(lambda x: x['name'].split('_')[0], axis=1)
            bench_df = utils.trim_benchmark_df(bench_df, col='true_name')
            bench_df = bench_df[~bench_df['sc_cst']]
            if args['--size']:
                sizes = df[args['--size']]
                sizes = normalize(sizes)
                kwargs['s'] = sizes
            else:
                kwargs['s'] = 10
            kwargs['c'] = colors.palette['red']
            kwargs['edgecolor'] = 'white'
            kwargs['linewidth'] = 0.05
            if args['--perres-x']:
                # bench_df[args['--xaxis']] = bench_df.apply(lambda x: x[args['--xaxis']] / x['chainA_size'], axis=1)
                bench_df[args['--xaxis']] = bench_df.apply(lambda x: x[args['--xaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)
            if args['--perres-y']:
                # bench_df[args['--yaxis']] = bench_df.apply(lambda x: x[args['--yaxis']] / x['chainA_size'], axis=1)
                bench_df[args['--yaxis']] = bench_df.apply(lambda x: x[args['--yaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)
            x = bench_df[args['--xaxis']]
            y = bench_df[args['--yaxis']]
            mpl_args = [x,y]
            points = ax.scatter(*mpl_args, marker='X', **kwargs)

        if args['--compare-final']:
            final_df = utils.safe_load(os.path.expanduser(
                '~/intelligent_design/helix_workspaces/design_models_final_combo_optimized/final.pkl'
            ))
            final_df = final_df[final_df['minimized']]
            print(final_df.columns)
            if args['--perres-x']:
                # final_df[args['--xaxis']] = final_df.apply(lambda x: x[args['--xaxis']] / x['chainA_size'], axis=1)
                final_df[args['--xaxis']] = final_df.apply(
                    lambda x: x[args['--xaxis']] / len(x['interfac_residues']), axis=1)
            if args['--perres-y']:
                # final_df[args['--yaxis']] = final_df.apply(lambda x: x[args['--yaxis']] / x['chainA_size'], axis=1)
                final_df[args['--yaxis']] = final_df.apply(
                    lambda x: x[args['--yaxis']] / len(x['interfac_residues']), axis=1)
            # if args['--target']:
            final_df = final_df[final_df['name'].isin(folder_dict[target_name])]
            x = final_df[args['--xaxis']]
            y = final_df[args['--yaxis']]
            mpl_args = [x, y]
            kwargs = {}
            if args['--size']:
                sizes = df[args['--size']]
                sizes = normalize(sizes)
                kwargs['s'] = sizes
            else:
                kwargs['s'] = 10
            # if hue:
            kwargs['edgecolor'] = 'white'
            kwargs['c'] = colors.palette['purple']
            kwargs['linewidth'] = 0.05
            points = ax.scatter(*mpl_args, picker=True, marker='d', **kwargs)
            # plt.colorbar(points, ax=ax)


        # ax.title.set_text(name_dict[target_name])
        # ax.set_title(name_dict[target_name], y=1.0, pad=-100, fontsize=8, loc='right')
        ax.text(0.97, 0.97, name_dict[target_name], horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
                fontsize=9)
        ax.tick_params(axis='x', labelsize=6, length=2)
        ax.tick_params(axis='y', labelsize=6, length=2)
        ax.tick_params(axis='both', which='major', pad=2)
    fig.text(0.5, 0.04, parse_metric(args['--xaxis']), fontsize=8, ha='center', va='top')
    fig.text(0.06, 0.5, parse_metric(args['--yaxis']), fontsize=8, ha='left', va='center', rotation='vertical')
    plt.subplots_adjust(
        left=0.125,
        right=0.9,
        bottom=0.1,
        top=0.9,
        wspace=0.2,
        hspace=0.2
    )

    # plt.tight_layout()
    plt.show()


def parse_dataframe(workspace, args):
    '''Parse dataframes with options'''
    print(workspace.target_path)
    df = workspace.get_scores(reread=args['--force-reread'])
    if 'chainA_size' not in df.columns:
        df = get_sizes(df)
        df.to_pickle(workspace.final_scorefile_path)
    if args['--filter']:
        print(' # designs before filtering:')
        print(df.shape[0])
        from helix.utils.utils import parse_filter
        import yaml
        with open(args['--filter'], 'r') as f:
            filter = yaml.load(f.read())
            filtered_df = parse_filter(filter, df)
        if not args['--keep-original']:
            df = filtered_df
            df['passed_filters'] = True
        else:
            df['passed_filters'] = False
            for idx in filtered_df.index:
                df.loc[idx, 'passed_filters'] = True
        print('# designs after filtering:')
        print(df.shape[0])
    if args['--model']:
        df['model_number'] = df.apply(lambda x: x['superimposed_file'].split('/')[-1].split('.')[0].split('_')[1], axis=1)
        df = df[df['model_number'] == args['--model']]
    df = df.reset_index()

    return df


def get_sizes(df):
    import prody
    sizes = {}
    size_col = []
    for idx, row in df.iterrows():
        model_num = os.path.basename(row['superimposed_file']).split('.')[0].split('_')[1]
        if model_num in sizes:
            size_col.append(sizes[model_num])
            continue
        atoms = prody.parsePDB(row['superimposed_file'], chain='A')
        sizes[model_num] = len(atoms.select('name CA').getIndices())
        size_col.append(sizes[model_num])
    df['chainA_size'] = pd.Series(size_col)
    return df

def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    if args['--target']:
        match_workspace = ws.MatchWorkspace(workspace.root_dir, args['--target'])
        workspaces = [match_workspace]
    else:
        workspaces = []
        targets = workspace.all_match_workspaces
        for target in targets:
            print(f'LOADING TARGET: {target}')
            if os.path.basename(target) == '3DI3':
                continue
            workspaces.append(ws.MatchWorkspace(workspace.root_dir, target))
    dfs = []
    for match_workspace in workspaces:
        df = parse_dataframe(match_workspace, args)

        if args['--suffix']:
            df = df[df['suffix'] == args['--suffix']]
        # print(df)
        # print(df.index)
        # print(df[df.total_score == df.total_score.max()])

        original_size = df.shape[0]
        df_dropped = df[~df['worst_9mer'].notna()]
        # print(df_dropped)
        # if df_dropped.shape[0] > 0:
            # print(df_dropped.design_file.iloc[0])
            # print(df_dropped.design_file.iloc[1])
        for col in [args['--xaxis'], args['--yaxis']]:
            df = df[df[col].notna()]
        no_dropped = original_size - df.shape[0]
        print(f'DF was originally {original_size} rows long')
        print(f'Dropped {no_dropped} rows due to missing values.')

        dfs.append(df)

    if args['<plot_type>'] == 'scatter':
        scatterplot(dfs, workspace, args)
    if args['<plot_type>'] == 'violin':
        violin_plot(dfs, workspace, args)


if __name__=='__main__':
    main()