'''
Run rifgen. Can split into multiple tasks if running on a job
distributor.

Usage:
    rifgen.py <folder> [options]

Options:
    --sge  Running on the cluster?
    --tasks=NUM, -n  How many tasks to split it into?
'''
import os, sys, docopt
import glob
import math
import subprocess
from subprocess import Popen, PIPE


def main():
    args = docopt.docopt(__doc__)
    print(args)
    folder = os.path.abspath(args['<folder>'])

    total_jobs = len(glob.glob(folder + '/*/'))
    #if '--tasks' in args:
    if args['--tasks']:
        num_tasks = int(args['--tasks'])
    else:
        num_tasks = 1

    #if '--sge' in args:
    if args['--sge']:
        task = os.environ['SGE_TASK_ID'] - 1
    else:
        task = 0
    
    start_job = task * math.ceil((total_jobs / num_tasks))
    stop_job = start_job + math.ceil(total_jobs / num_tasks) - 1

    folders = sorted(glob.glob(folder + '/*/'))

    rifgen = os.environ['RIFGEN']
    for fold in folders[start_job:stop_job+1]:
        print(fold)
        os.chdir(fold)
        flags = os.path.join(fold, 'flags')
        process = Popen([rifgen, '@', flags], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        exit_code = process.wait()
        if exit_code != 0:
            print(stdout)
            print(stderr)
        output = subprocess.check_output([rifgen, '@', flags])
        print(output)


if __name__=='__main__':
    main()
