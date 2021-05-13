import pandas as pd
import glob
import sys
import os
import pickle5 as pickle

out = pd.DataFrame()
for path in sorted(glob.glob(sys.argv[1] + '/*.pkl')):
    print(path)
    # df = pd.read_pickle(path)
    with open(path, 'rb') as f:
        df = pickle.load(f)
    out = out.append(df, ignore_index=True)

out.to_pickle(sys.argv[1] + '/final.pkl')
