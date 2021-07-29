'''
Quick command that combines pickled dataframes from different tasks into a
single, "final" dataframe.

Usage:
    helix combine_dataframes <folder> [options]

Options:
    --out=FILE, -o  Where to save the final dataframe
'''


import pandas as pd
import docopt
import glob
import sys
import os


args = docopt.docopt(__doc__)
out = pd.DataFrame()
for path in sorted(glob.glob(args['<folder>'] + '/*.pkl')):
    print('Reading {}'.format(path))
    df = pd.read_pickle(path)
    out = out.append(df, ignore_index=True)

if args['--out']:
    out.to_pickle(args['--out'])
else:
    out.to_pickle(sys.argv[1] + '/final.pkl')
