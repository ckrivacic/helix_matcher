"""
Usage:
    view_candidates.py <candidates_dataframe> [options]

Options:
    --final_df=PKL, -f  Start where you left off in the final dataframe.
    --out=PATH  Where to save the final candidates  
    [default: interface_finder/final_candidates.pkl]
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

    if args['--final_df']:
        final = pd.read_pickle(args['--final_df'])
        max_idx = max(final.index)
        df = df.loc[max_idx + 1:]
    script = os.path.abspath("launch_pymol_candidates.sho")
    yes = ['y', 'yes', 'oui', 'si', 'yeppers', 'yessir', 'heckyeah']
    no = ['no', 'n', 'non', 'noway', 'fohgetaboudit']
    out = args['--out']
    groups = df.groupby(['name', 'target', 'chain'])
    total_groups = groups.ngroups
    i = 0
    for name, group in groups:
        i += 1
        print("Viewing group {} of {}".format(i, total_groups))
        selstr = "sele helices, {} and (".format(name[0])
        for idx, row in group.iterrows():
            selstr += "(resi {}-{} and chain {}) or ".format(
                    min(row['pdb_resis']), max(row['pdb_resis']),
                    row['chain']
                    )
        selstr = selstr[:-4] + ')'
        print(selstr)
        cmd = script, name[0], selstr, row['chain'], row['target']
        subprocess.call(cmd)
        save_example = input("Save this benchmark candidate? (y/N)") 
        if save_example in yes:
            if os.path.exists(out):
                final_df = pd.read_pickle(out)
                final_df = pd.concat([final_df, group])
                final_df.to_pickle(out)
            else:
                group.to_pickle(out)


if __name__=='__main__':
    main()
