'''
Run rifgen. Can split into multiple tasks if running on a job
distributor.

Usage:
    rifgen.py <folder> [Options]

Options:
    --sge  Running on the cluster?
    --tasks=NUM, -n  How many tasks to split it into?
'''
import os, sys, docopt
import glob
import math
from subprocess import Popen


def main():
    args = docopt.docopt(__doc__)
    folder = args['<folder>']

    total_jobs = len(glob.glob(folder + '/*/'))
    if args['--tasks']:
        num_tasks = int(args['--tasks'])
    else:
        num_tasks = total_jobs

    if args['--sge']:
        task = os.environ['SGE_TASK_ID'] - 1
    else:
        task = 0
    
    start_job = task * math.ceil((total_jobs / num_tasks))
    stop_job = start_job + math.ceil(total_jobs / num_tasks) - 1

    folders = sorted(glob.glob(folder + '/*/'))

    rifgen = os.environ['RIFGEN']
    for fold in folders[start_job:stop_job+1]:
        flags = os.path.join(fold, 'flags')
        Popen([rifgen, '@{}'.format(flags)])
