'''
Analyze different design methods.

Usage:
    design_analysis.py <plot_type>
'''
import docopt
from helix.utils import utils


def calc_identity(row):
    '''Calculate sequence identity of a helix compared to the nearest
    benchmark helix'''



def plot_sequence_recovery(df, args):
    '''Make a plot for sequence recovery'''
    df = df[df['best_rmsd'] < 1.0]
    df['helix_identity'] = df.apply()
    groups = df.groupby('focus_dir')


def main():
    args = docopt.docopt(__doc__)
    plot_type = args['<plot_type>']
    df = utils.safe_load('rifdock_outputs/combined_benchmark_rev/final.pkl')
    
    if plot_type == 'seq':
        plot_sequence_recovery(df, args)


if __name__=='__main__':
    main()
