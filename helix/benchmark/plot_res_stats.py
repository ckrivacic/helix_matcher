'''
Plot interface residue statistics such as # of helical residues vs. other secstructs.

Usage:
    plot_res_stats.py [options]

Options:
'''

import docopt
from helix.utils import utils
import os
import seaborn as sns
from numpy import median
import pandas as pd
from helix.utils.colors import palette
from matplotlib import pyplot as plt
import matplotlib as mpl


def make_barplot(df, x, y, hue=None):
    colors = [palette['lightgray'], palette['darkgray'], 'black']
    order = ['surface', 'boundary', 'buried']
    hue_order = ['H', 'E', 'L']

    fig, ax = plt.subplots(figsize=(4,3), dpi=300)
    sns_ax = sns.barplot(data=df, x=x, y=y,
                         hue=hue, order=order, hue_order=hue_order, ax=ax, saturation=1, palette=colors,
                         estimator=median,
                        ci=95)
    for patch in sns_ax.patches:
        patch.set_edgecolor('black')

    labels = ['Surface', 'Boundary', 'Buried']
    ax.set_xticklabels(labels, ha='right')
    plt.xticks(rotation=70)
    plt.ylabel(y)
    plt.xlabel(None)

    plt.show()


def main():
    mpl.use('tkagg')
    dfpath = os.path.join('residue_scores_min/final.pkl')
    df = utils.safe_load(dfpath)
    df = df.rename(columns={'burial': 'Burial', 'secstruct': 'Secondary structure', 'total_crosschain': 'Median intermolecular REU'})
    print(df.columns)
    x = 'Burial'
    hue = 'Secondary structure'
    y = 'Median intermolecular REU'

    make_barplot(df, x, y, hue=hue)


    rows = []
    for xval in df[x].unique():
        df_x = df[df[x] == xval]
        total_size = df_x.shape[0]
        for h in df_x[hue].unique():
            df_hue = df_x[df_x[hue] == h]
            fraction = df_hue.shape[0] / total_size
            row = {
                'Burial': xval,
                'Secondary structure': h,
                'Fraction': fraction,
            }
            rows.append(row)

    df = pd.DataFrame(rows)
    make_barplot(df, 'Burial', 'Fraction', 'Secondary structure')

if __name__=='__main__':
    main()