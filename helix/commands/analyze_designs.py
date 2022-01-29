'''
Analyze different design methods.

Usage:
    helix analyze_designs <workspace> <plot_type> [options]

Options:
    --yaxis=COLUMN, -y  Which column to plot on the y-axis  
        [default: helix_percent_identity]

    --xaxis=COLUMN, -x  Which column to plot on the x-axis  
        [default: interface_score_y]

    --rmsd-cutoff=FLOAT, -c  Eliminate helices greater than this RMSD cutoff  
        [default: 1.0]

    --pose-score-cutoff=FLOAT, -p  Elimiate poses with a total score
        greater than this float  [default: 0]

    --focus-dir=STR, -f  Only plot for a specific focus dir or list of focus
        dirs (ex. 1b33_K,1b33_K_specialrot)

    --size=COLUMN, -s  Size scatter plot points by this column

    --hue=COLUMN, -h  Color by a column  [default: focus_dir]
'''
import docopt
import os
import pandas as pd
from helix.utils import utils
import helix.workspace as ws
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt


def calc_identity(row):
    '''Calculate sequence identity of a helix compared to the nearest
    benchmark helix'''
    seq_helix = row['helix_seq']
    seq_benchmark = row['benchmark_seq']
    assert len(seq_helix) == len(seq_benchmark)
    identity = 0
    for i, res in enumerate(seq_helix):
        if res == seq_benchmark[i]:
            identity += 1

    return 100 * (identity / len(seq_helix))


def plot_sequence_recovery(df, args):
    '''Make a plot for sequence recovery'''
    # for col in df.columns:
        # printcol = col
        # print(df.sort_values(by=printcol)[printcol])
    order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            '1b33_K_nativelike', '1b33_K_specialrot']
    sns.swarmplot(data=df, x='focus_dir', y=args['--yaxis'],
            order=order, color='.5', alpha=0.5)
    sns.violinplot(data=df, x='focus_dir', y=args['--yaxis'],
            order=order)
    plt.show()


def scatterplot(df, args):
    if args['--size']:
        size = args['--size']
        sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=args['--hue'], size=size, sizes=(50,150))
    else:
        sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=args['--hue'])
    plt.show()


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    plot_type = args['<plot_type>']
    dfpath = os.path.join(workspace.rifdock_outdir,
            'combined_benchmark_rev', 'final.pkl')
    df = utils.safe_load(dfpath)
    if args['--focus-dir']:
        focusdirs = args['--focus-dir'].split(',')
        df = df[df['focus_dir'].isin(focusdirs)]
    print('Possible values for x- and y-axes are:')
    for col in df.columns:
        print(col)
    if args['--yaxis'] == 'helix_percent_identity':
        df['helix_percent_identity'] = df.apply(calc_identity, axis=1)
        df.to_pickle('temp.pkl')
        df = pd.read_pickle('temp.pkl')
    if not float(args['--rmsd-cutoff']) == 0:
        df = df[df['best_rmsd'] < float(args['--rmsd-cutoff'])]
    if not args['--pose-score-cutoff'] == 'False':
        df = df[df['pose_score'] < float(args['--pose-score-cutoff'])]
    
    if plot_type == 'seq':
        plot_sequence_recovery(df, args)
    if plot_type == 'scatter':
        scatterplot(df, args) 


if __name__=='__main__':
    main()
