'''
Plot forward folding results from a Roseasy workspace

Usage:
    helix plot_ff_roseasy <roseasy_workspace> [options]

Options:
    --plot-type=STR  Which type of plot to produce, options: scatter,  [default: scatter]
    --xaxis=COLUMN, -x  Plot this column on the x-axis  [default: ca_rmsd]
    --yaxis=COLUMN, -y  Plot this column on the y-axis  [default: holes_monomer]
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
import glob, os
import yaml
import matplotlib as mpl
from roseasy import pipeline


def get_scores(folder, reread=False):
    # This will go faster if you run helix combine on the "scores"
    # folder
    # Reread param makes you open individual datafiles regardless, combining any that are
    # not present in the final dataframe.
    final_scorefile_path = os.path.join(folder, 'final.pkl')
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


def parse_dataframe(df, workspace, args):
    '''Parse dataframes with options'''
    final_scorefile_path = os.path.join(workspace.focus_dir, 'final.pkl')

    def relative_path(row):
        path = row['design']
        pathlist = path.split('/')
        root_idx = pathlist.index(os.path.basename(workspace.root_dir)) + 1
        relpath = os.path.join(workspace.root_dir, *pathlist[root_idx:])
        return relpath

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
    for dir in focus_dirs:
        focus_ws = pipeline.workspace_from_dir(dir)
        score_dir = os.path.join(dir, 'analysis')
        df = get_scores(score_dir)
        df = parse_dataframe(df, focus_ws, args)

        original_size = df.shape[0]
        for col in [args['--xaxis'], args['--yaxis']]:
            df = df[df[col].notna()]
        no_dropped = original_size - df.shape[0]
        print(f'DF was originally {original_size} rows long')
        print(f'Dropped {no_dropped} rows due to missing values.')

        dfs.append(df)

    if args['--plot-type'] == 'scatter':
        scatterplot(dfs, workspace, args)