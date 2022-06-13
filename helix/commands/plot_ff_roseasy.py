'''
Plot forward folding results from a Roseasy workspace

Usage:
    helix plot_ff_roseasy <roseasy_workspace> [options]

Options:
    --plot-type=STR  Which type of plot to produce, options: scatter,  [default: scatter]
    --xaxis=COLUMN, -x  Plot this column on the x-axis  [default: ca_rmsd]
    --yaxis=COLUMN, -y  Plot this column on the y-axis  [default: npsa_monomer]
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
import docopt
from helix.commands.plot_designs import scatterplot
from helix.utils.utils import parse_filter
from helix.commands.plot_designs import get_sizes
from helix.utils.utils import safe_load
import pandas as pd
import numpy as np
from pyrosetta import *
import seaborn as sns
import glob, os
import yaml
import matplotlib as mpl
import matplotlib.pyplot as plt
from roseasy import pipeline


def get_scores(folder, reread=False):
    # This will go faster if you run helix combine on the "scores"
    # folder
    # Reread param makes you open individual datafiles regardless, combining any that are
    # not present in the final dataframe.
    final_scorefile_path = os.path.join(folder, 'metrics.pkl')
    if os.path.exists(final_scorefile_path) and not reread:
        # if reread:
        #     old_df = safe_open_dataframe(self.final_scorefile_path)
        # else:
        return safe_load(final_scorefile_path)

    dataframes = []
    all_scorefiles = glob.glob(os.path.join(folder, '*.pkl'))
    for f in all_scorefiles:
        df = safe_load(f)
        dataframes.append(df)
    if len(dataframes) > 0:
        df = pd.concat(dataframes, ignore_index=True)
    else:
        return pd.DataFrame()
    # if reread:
    #     df = pd.concat([df, old_df]).drop_duplicates('design_file').reset_index(drop=True)
    df.to_pickle(final_scorefile_path)
    return df


def scatterplot(df, workspace, args, use_matplotlib=True):
    groups = df.groupby('design_name')
    groupkeys = list(groups.groups.keys())
    # for idx, ax in enumerate(axs.reshape(-1)):
    sfxn = create_score_function('ref2015')
    for name, df in groups:
        original_pose = df.iloc[0]['design_relpath']
        original_pose = os.path.join(workspace.root_dir, original_pose)
        original_pose = pose_from_file(original_pose)
        score = sfxn(original_pose)
        new_row = {
            'total_score': score,
            'ca_rmsd': 0,
            'target': 'input',
        }
        orig_df = pd.DataFrame([new_row])
        df = pd.concat([df, orig_df], ignore_index=True)
        print(df)
        # df = groups.groups[groupkeys[idx]]
        # if idx > len(groups)-1:
        #     continue
        print('DF COLS')
        print(df.columns)
        # if df.shape[0] > 0:
        #     target_name = df.iloc[0]['design_file'].split('/')[1]
        # else:
        #     target_name = 'N/A'
        if args['--perres-x']:
            # df[args['--xaxis']] = df.apply(lambda x: x[args['--xaxis']] / x['chainA_size'], axis=1)
            df[args['--xaxis']] = df.apply(lambda x: x[args['--xaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)
        if args['--perres-y']:
            # df[args['--yaxis']] = df.apply(lambda x: x[args['--yaxis']] / x['chainA_size'], axis=1)
            df[args['--yaxis']] = df.apply(lambda x: x[args['--yaxis']] / len(x['interfac_residues']) if len(x['interfac_residues']) > 0 else 0, axis=1)

        # if use_matplotlib:
        #     # Seaborn seems to reorder points, making the on_pick stuff useless. Back to MPL.
        #     x = df[args['--xaxis']]
        #     y = df[args['--yaxis']]
        #     mpl_args = [x, y]
        #     kwargs = {}
        #     if args['--size']:
        #         sizes = df[args['--size']]
        #         sizes = normalize(sizes)
        #         kwargs['s'] = sizes
        #     else:
        #         kwargs['s'] = 1
        #     if hue:
        #         kwargs['c'] = df[hue]
        #     else:
        #         kwargs['c'] = colors.palette['teal']
        #     kwargs['alpha'] = 0.5
        #     kwargs['edgecolor'] = 'white'
        #     kwargs['linewidth'] = 0.05
        #     points = ax.scatter(*mpl_args, picker=True, **kwargs)
        #     if args['--hue']:
        #         plt.colorbar(points, ax=ax)
        #     click = plotting.ClickablePlot(points, df, args, workspace)
        # else:
        ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                             alpha=1.0, hue='target')
        if args['--yaxis'] == 'total_score':
            max_y_df = df[df['total_score'] < 0]
            max_y = max(max_y_df['total_score'])
            ax.set_ylim(min(df[args['--yaxis']]), max_y)

        ax.title.set_text("{}_{}".format(df.iloc[0]['target'], df.iloc[0]['design_name']))
        # ax.set_title(name_dict[target_name], y=1.0, pad=-100, fontsize=8, loc='right')
        # ax.text(0.97, 0.97, name_dict[target_name], horizontalalignment='right', verticalalignment='top', transform=ax.transAxes,
        #         fontsize=9)
        # ax.tick_params(axis='x', labelsize=6, length=2)
        # ax.tick_params(axis='y', labelsize=6, length=2)
        # ax.tick_params(axis='both', which='major', pad=2)
    # fig.text(0.5, 0.04, parse_metric(args['--xaxis']), fontsize=8, ha='center', va='top')
    # fig.text(0.06, 0.5, parse_metric(args['--yaxis']), fontsize=8, ha='left', va='center', rotation='vertical')
    # plt.subplots_adjust(
    #     left=0.125,
    #     right=0.9,
    #     bottom=0.1,
    #     top=0.9,
    #     wspace=0.2,
    #     hspace=0.2
    # )

    # plt.tight_layout()
        plt.show()


def parse_dataframe(df, workspace, args):
    '''Parse dataframes with options'''
    final_scorefile_path = os.path.join(workspace.focus_dir, 'final.pkl')

    def relative_path(row):
        path = row['design']
        pathlist = path.split('/')
        root_idx = pathlist.index(os.path.basename(workspace.root_dir)) + 1
        relpath = os.path.join(workspace.root_dir, *pathlist[root_idx:])
        return relpath

    print(df)
    df['design_realpath'] = df.apply(relative_path, axis=1)

    if not 'model_number' in df.columns:
        df['model_number'] = df.apply(lambda x: x['design_name'].split('/')[-1].split('.')[0].split('_')[1], axis=1)
    if args['--model']:
        df = df[df['model_number'] == args['--model']]

    if 'chainA_size' not in df.columns:
        df = get_sizes(df, name_col='design_realpath')
        df.to_pickle(final_scorefile_path)
    if args['--filter']:
        print(' # designs before filtering:')
        print(df.shape[0])
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
    df = df.reset_index()

    return df


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    workspace = pipeline.workspace_from_dir(args['<roseasy_workspace>'])

    focus_dirs = glob.glob(os.path.join(workspace.root_dir, '*_validated_designs'))
    dfs = []
    init()
    for dir in focus_dirs:
        print(f"Loading from dir {dir}")
        focus_ws = pipeline.workspace_from_dir(dir)
        score_dirs = glob.glob(os.path.join(focus_ws.output_dir, '*/*/'))
        # score_dirs = glob.glob(os.path.join(focus_ws.focus_dir, 'analysis/'))
        for score_dir in score_dirs:
            df = get_scores(score_dir)
            if df.empty:
                continue
            df = parse_dataframe(df, focus_ws, args)

            original_size = df.shape[0]
            for col in [args['--xaxis'], args['--yaxis']]:
                df = df[df[col].notna()]
            no_dropped = original_size - df.shape[0]
            print(f'DF was originally {original_size} rows long')
            print(f'Dropped {no_dropped} rows due to missing values.')

            print('Dropping positive scores')

            # df = df[df['total_score'] < 0]

            dfs.append(df)
            scatterplot(df, workspace, args)

    if args['--plot-type'] == 'scatter':
        for df in dfs:
            scatterplot(df, workspace, args)