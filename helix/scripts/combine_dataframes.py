import pandas as pd
import glob
import sys
import os

out = pd.DataFrame()
for path in sorted(glob.glob(sys.argv[1] + '/*.pkl')):
    print(path)
    df = pd.read_pickle(path)
    out = out.append(df, ignore_index=True)

out.to_pickle(sys.argv[1] + '/final.pkl')
