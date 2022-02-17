import pandas as pd
import os
import sys
import glob


'''
combine scoring dataframes with benchmark dataframes
'''


def update_benchmark(scoredf, benchmark_df, reverse=True):
    '''Update the benchmark dataframe with scoring info.'''
    if reverse:
        return pd.merge(benchmark_df, scoredf, on=['patchman_file'],
                validate='1:1')
    else:
        scoredf.rename(columns={'patchman_file':'len_14_path'},
                inplace=True)
        return pd.merge(benchmark_df, scoredf, on=['len_14_path'],
                validate='m:m')
    # print(benchmark_df)
    # for idx, row in scoredf.iterrows():
        # benchmark_row = benchmark_df[benchmark_df==row['name']]


# scoretype = 'pose_score'
scoretype = 'interface_score'

benchmark_dfpath = 'benchmark_results/final.pkl'
benchmark_reverse_dfpath = 'benchmark_results_reverse/final.pkl'
# dfpath = benchmark_reverse_dfpath
dfpath = benchmark_reverse_dfpath
benchmark_df = pd.read_pickle(dfpath)
def func(x):
    return os.path.basename(x)
# print('Appending name')
# benchmark_df['name'] = benchmark_df.apply(lambda x:
        # func(x['patchman_file']), axis=1)
# print(benchmark_df['name'])

# length = 'len_14'
# globstr = 'rifdock_outputs/*/patch_*/{}/docked_full/*.pdb.gz'.format(length)

# for filename in glob.glob(globstr):
    # if benchmark_dfpath = ''
    # row = benchmark_df[benchmark_df[]]

scoreglob = sorted(glob.glob('*/scores/task_*.pkl'))
out_dfs = []
# Task 1-indexed
if len(sys.argv) > 1:
    task = int(sys.argv[1]) - 1
else:
    task = int(os.environ['SGE_TASK_ID']) - 1

total_tasks = len(scoreglob)
interval = 500
start = task * interval
for scorefile in scoreglob[start:start + interval]:
    print('Updating benchmark dataframe for scores in {}'.format(
        scorefile
        ))
    scoredf = pd.read_pickle(scorefile)
    out_dfs.append(update_benchmark(scoredf, benchmark_df,
        reverse=True))

final = pd.concat(out_dfs, ignore_index=True)
final.to_pickle('combined_benchmark_rev/task_{}.pkl'.format(task))
