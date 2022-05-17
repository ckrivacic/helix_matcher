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
    --split-by-target  Plot all target on one figure, but with a separate plot for each target
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


def scatterplot(df, workspace, args):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
    # '1b33_K_nativelike', '1b33_K_specialrot']
    if not args['--xaxis'] == 'protocol':
        order = None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    print(df[args['--yaxis']].sort_values())
    if args['--size']:
        size = args['--size']
        ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                             hue=hue, size=size, sizes=(50,300), picker=True)
    else:
        ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                             hue=hue, picker=True)
    click = plotting.ClickablePlot(ax, df, args, workspace)

    plt.show()


def parse_dataframe(workspace, args):
    '''Parse dataframes with options'''
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
            workspaces.append(ws.MatchWorkspace(workspace.root_dir, target))
    for match_workspace in workspaces:
        df = parse_dataframe(match_workspace, args)
        print(df)

    if args['<plot_type>'] == 'scatter':
        scatterplot(df, workspace, args)


if __name__=='__main__':
    main()