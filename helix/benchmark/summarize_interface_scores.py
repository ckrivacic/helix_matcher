'''
Take residue data and collate it into something usable as a score.
Optionaly plot distributions.

Usage:
    summarize_interface_scores.py [options]

Options:
    --dataframe=PATH, -d  Path to the dataframe containing all interface
        residue information.  [default: residue_scores/final.pkl]
    --by=COMMA_SEP_LIST  Which columns the data should be grouped by
        prior to summarization. For multiple columns, input a
        comma-separated list (no spaces), i.e. restype,burial  [default: restype]
    --plot=DATATYPE, -p  Plot data of a particular type (column)
        [default: total_energy]
    --plot-by=COMMA_SEP_LIST  Which column should go into separate
        plots of plotting  [default: restype]
    --cat=CATEGORY  Color by this category  [default: burial]
    --crosschain  Plot crosschain summary
    --fa  Plot fullatom summary
'''

import pickle5 as pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import docopt
import numpy as np
from scipy import stats
from helix.utils import utils


def summarize_data(df, group_cols):
    '''Get mean, median, and SD for each group'''
    scoretypes = [
            'fa_atr_cc', 'fa_rep_cc', 'fa_sol_cc', 'fa_elec_cc',
            'hbond_bb_sc_cc', 'hbond_sc_cc', 'total_crosschain',
            'fa_atr_tot', 'fa_rep_tot', 'fa_sol_tot', 'fa_elec_tot',
            'hbond_bb_sc_tot', 'hbond_sc_tot', 'pro_close_tot',
            'fa_intra_rep_tot', 'dslf_fa13_tot', 'rama_tot',
            'omega_tot', 'fa_dun_tot', 'p_aa_pp_tot', 'ref_tot',
            'hbond_sr_bb_tot', 'total_energy',
            ]
    if 'contacts' not in group_cols:
        scoretypes.append('contacts')
    groups = df.groupby(group_cols)
    outrows = []
    for name, group in groups:
        row0 = group.iloc[0]
        row = {}
        for col in group_cols:
            row[col] = row0[col]
        for st in scoretypes:
            subrow = row.copy()
            subrow['scoretype'] = st
            subrow['mean'] = group[st].mean()
            subrow['median'] = group[st].median()
            subrow['sd'] = group[st].std()
            outrows.append(subrow)

    return pd.DataFrame(outrows)


def plot_dist(df, args):
    '''
    Plot distributions (not summarized data)
    '''
    plot_by = args['--plot-by'].split(',')
    # df = df[(np.abs(stats.zscore(df[args['--plot']])) < 3)]
    # df = df[(np.abs(stats.zscore(df[args['--plot']])) < 3)]
    # df = df[df[args['--plot']] < 10]
    groups = df.groupby(plot_by)
    if plot_by[0]=='restype':
        fig, axes = plt.subplots(4, 5, sharey='all',)
        fig.set_figheight(15)
        fig.set_figwidth(15)
        # plt.ylim(0, 1)
    else:
        fig, axes = plt.subplots(groups.ngroups)
    cat = args['--cat']
    alpha = 1 / len(set(df[cat]))
    datrange = max(df[args['--plot']]) - min(df[args['--plot']])
    for (name, group), ax in zip(groups, axes.flatten()):
        ax.set_title(name)
        for subname, subgroup in group.groupby(cat):
            sns.distplot(subgroup[args['--plot']], kde_kws={'linewidth': 1},
                    ax=ax, bins=int(datrange / 0.5), label=subname)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_summarized(df, args):
    '''
    plot summarized data
    '''
    plot_by = args['--plot-by'].split(',')
    # df = df[df['scoretype']==args['--plot']]
    # df = df[(np.abs(stats.zscore(df[args['--plot']])) < 3)]
    # df = df[(np.abs(stats.zscore(df[args['--plot']])) < 3)]
    df = df[df[args['--plot']] < 10]
    sns.barplot(data=df, x=args['--plot-by'], y=args['--plot'],
            hue=args['--cat'],)
    plt.show()


def plot_all_scores(df, args):
    cc_scoretypes = [
            'fa_atr_cc', 'fa_rep_cc', 'fa_sol_cc', 'fa_elec_cc',
            'hbond_bb_sc_cc', 'hbond_sc_cc', 'total_crosschain',]
    fa_scoretypes = [
            'fa_atr_tot', 'fa_rep_tot', 'fa_sol_tot', 'fa_elec_tot',
            'hbond_bb_sc_tot', 'hbond_sc_tot', 'pro_close_tot',
            'fa_intra_rep_tot', 'dslf_fa13_tot', 'rama_tot',
            'omega_tot', 'fa_dun_tot', 'p_aa_pp_tot', 'ref_tot',
            'hbond_sr_bb_tot', 'total_energy',
            ]
    if args['--crosschain']:
        df = df[df['scoretype'].isin(cc_scoretypes)]
    if args['--fa']:
        df = df[df['scoretype'].isin(fa_scoretypes)]
    # for scoretype in scoretypes:
        # df = df[(np.abs(stats.zscore(df[scoretype])) < 3)]
    sns.barplot(data=df, x='restype', y='mean', hue='scoretype')
    plt.show()


def main():
    canonical_aas = ['GLU', 'CYS', 'HIS', 'GLY', 'GLN', 'ALA', 'ASN', 'PRO', 'MET', 'TYR', 'PHE', 'ASP', 'TRP', 'LEU', 'SER', 'THR', 'LYS', 'ILE', 'VAL', 'ARG']
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    print('Loading dataframe...')
    df = utils.safe_load(args['--dataframe']) 
    print('Dataframe loaded. Triming noncanonical AAs.')
    df = df[df['restype'].isin(canonical_aas)]
    print('Noncanonical AAs trimmed. Grouping dataframe.')
    group_cols = args['--by'].split(',')
    print('Dataframe grouped. Summarizing dataframe.')
    summarized = summarize_data(df, group_cols)
    print('Dataframe summarized.')
    print(summarized)
    # plot_summarized(df, args)
    plot_dist(df, args)
    # plot_all_scores(summarized, args)


if __name__=='__main__':
    main()
