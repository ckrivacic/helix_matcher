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
    --all-targets, -a  Plot all target on one figure, but with a separate plot for each target
    --force-reread, -f  Even if final.pkl exists in the design directory, reread all individual score files.
    --filter=FILEPATH  Filter the dataframe
    --keep-original  Keep original dataframe when filtering, plot both
'''
from helix.utils import utils
import helix.workspace as ws
from helix.utils import plotting
import docopt
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math


def normalize(array, min_val=3, max_val=50):
    array = array - min_val
    x = max(array) / max_val
    array = array / x
    return array


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
        nrows = 2
        ncols = math.ceil(len(dfs) / 2)
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols)
    axs = np.array(axs)
    # if args['--size']:
    #     size = args['--size']
    #     ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
    #                          hue=hue, size=size, sizes=(50,300), picker=True)
    for idx, ax in enumerate(axs.reshape(-1)):
        if idx > len(dfs)-1:
            continue
        df = dfs[idx]
        target_name = df.iloc[0]['design_file'].split('/')[1]
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
                kwargs['s'] = 20
            if hue:
                kwargs['c'] = df[hue]
            points = ax.scatter(*mpl_args, picker=True, **kwargs)
            plt.colorbar(points, ax=ax)
            ax.title.set_text(target_name)
            click = plotting.ClickablePlot(points, df, args, workspace)
        else:
            ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                                 hue=hue, picker=True)

    plt.show()


def parse_dataframe(workspace, args):
    '''Parse dataframes with options'''
    print(workspace.target_path)
    df = workspace.get_scores(reread=args['--force-reread'])
    if args['--filter']:
        from helix.commands.filter import parse_filter
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
    df = df.reset_index()

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
        print(df_dropped)
        print(df_dropped.design_file.iloc[0])
        print(df_dropped.design_file.iloc[1])
        for col in [args['--xaxis'], args['--yaxis']]:
            df = df[df[col].notna()]
        no_dropped = original_size - df.shape[0]
        print(f'Dropped {no_dropped} rows due to missing values.')
        dfs.append(df)

    if args['<plot_type>'] == 'scatter':
        scatterplot(dfs, workspace, args)


if __name__=='__main__':
    main()