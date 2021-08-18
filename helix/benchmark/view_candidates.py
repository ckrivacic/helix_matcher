"""
Usage:
    view_candidates.py <candidates_dataframe>
"""


import pandas as pd
import pickle5 as pickle
import os
import docopt
import subprocess


def main():
    args = docopt.docopt(__doc__)
    dataframe_path = args['<candidates_dataframe>']
    try:
        df = pd.read_pickle(dataframe_path)
    except:
        with open(dataframe_path, 'rb') as f:
            df = pickle.load(f)

    script = os.path.abspath("launch_pymol_candidates.sho")
    for name, group in df.groupby(['name']):
        selstr = "sele helices, {} and (".format(name)
        for idx, row in group.iterrows():
            selstr += "(resi {}-{} and chain {}) or ".format(
                    min(row['pdb_resis']), max(row['pdb_resis']),
                    row['chain']
                    )
        selstr = selstr[:-4] + ')'
        print(selstr)
        cmd = script, name, selstr
        subprocess.call(cmd)


if __name__=='__main__':
    main()
