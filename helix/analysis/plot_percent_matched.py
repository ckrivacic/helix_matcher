'''
plot percent of helices and scaffolds lost at each step

Usage:
    plot_percent_matched.py <workspace> [options]

Options:
'''
from helix.utils import utils
import helix.workspace as ws
import docopt
import os
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from helix.utils.colors import palette
import matplotlib.pyplot as plt


def barplot(df, args):
    fig, ax = plt.subplots(figsize=(4,3), dpi=300)
    ax1 = sns.barplot(data=df, x='Protein', y='s_fraction_matched', palette=[palette['teal']], ax=ax, ci=None)
    ax2 = sns.barplot(data=df, x='Protein', y='s_fraction_scored', palette=[palette['blue']], ax=ax, ci=None)
    ax3 = sns.barplot(data=df, x='Protein', y='s_fraction_designed', palette=[palette['darkgray']], ax=ax, ci=None)
    delt = 0.135
    width = 0.4
    for axis in [ax1, ax2, ax3]:
        for i, bar in enumerate(axis.patches):
            bar.set_x(bar.get_x() - delt)
            bar.set_width(width)
            bar.set_edgecolor('black')

    # labels = df['Total scaffolds']
    # for rect, label in zip(ax1.patches, labels):
    #     height = rect.get_height()
    #     ax.annotate(label,
    #          (rect.get_x() - delt + 0.7, height + 0.01), ha='center', va='bottom',
    #                 fontsize=4
    #     )


    ax1 = sns.barplot(data=df, x='Protein', y='h_fraction_matched', palette=[palette['orange']], ax=ax, ci=None)
    ax2 = sns.barplot(data=df, x='Protein', y='h_fraction_scored', palette=[palette['red']], ax=ax, ci=None)
    ax3 = sns.barplot(data=df, x='Protein', y='h_fraction_designed', palette=[palette['darkgray']], ax=ax, ci=None)
    # ax.bar_label(ax1.containers[-3], fmt="Total:\n%", labels=df['total # helices'], fontsize=8)
    for axis in [ax1, ax2, ax3]:
        for i, bar in enumerate(axis.patches):
            bar.set_x(bar.get_x() + delt)
            bar.set_width(width)
            bar.set_edgecolor('black')

    # labels = df['total # helices']
    # for rect, label in zip(ax1.containers[-3], labels):
    #     height = rect.get_height()
    #     ax1.annotate(label,
    #                  (rect.get_x(), height + 0.04), ha='center', va='bottom',
    #                  fontsize=6
    #                  )

    ax.set_xticklabels(df['Protein'], ha='right')
    plt.xticks(rotation=70)
    plt.ylabel('Fraction scaffolds/helices')

    plt.show()


def calc_percent(row, col1, col2):
    return row[col1] / row[col2]

def main():
    args = docopt.docopt(__doc__)
    mpl.use('tkagg')
    workspace = ws.workspace_from_dir(args['<workspace>'])
    df = pd.read_csv(os.path.join(workspace.root_dir, 'helix_percent_data.csv'))
    df['h_fraction_matched'] = df.apply(calc_percent, args=('# helices matched', 'total # helices'), axis=1)
    df['h_fraction_scored'] = df.apply(calc_percent, args=('# helices passed basic filters', 'total # helices'), axis=1)
    df['h_fraction_designed'] = df.apply(calc_percent, args=('# designed helices', 'total # helices'), axis=1)
    df['s_fraction_matched'] = df.apply(calc_percent, args=('# scaffolds matched', 'Total scaffolds'), axis=1)
    df['s_fraction_scored'] = df.apply(calc_percent, args=('# scaffolds passed basic filters', 'Total scaffolds'), axis=1)
    df['s_fraction_designed'] = df.apply(calc_percent, args=('# scaffolds designed', 'Total scaffolds'), axis=1)

    barplot(df, args)



if __name__=='__main__':
    main()