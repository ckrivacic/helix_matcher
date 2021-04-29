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

    total_jobs = len(glob.glob(folder + '/patch_*'))
    print('TOTAL JOBS: {}'.format(total_jobs))
    #if '--tasks' in args:
    if args['--tasks']:
        num_tasks = int(args['--tasks'])
    else:
        num_tasks = total_jobs

    #if '--sge' in args:
    if args['--sge']:
        task = int(os.environ['SGE_TASK_ID']) - 1
    else:
        task = 0

    print('TASK: {}'.format(task))
    
    start_job = task * math.ceil((total_jobs / num_tasks))
    stop_job = start_job + math.ceil(total_jobs / num_tasks)
    print('START JOB: {}'.format(start_job))
    print('STOP JOB: {}'.format(stop_job))

    folders = sorted(glob.glob(folder + '/patch_*'))

    rifgen = os.path.join(folder, 'rifgen')
    for fold in folders[start_job:stop_job]:
        print(fold)
        os.chdir(fold)
        flags = os.path.join(fold, 'flags')
        myenv = os.environ.copy()
        myenv['LD_LIBRARY_PATH'] = '/wynton/home/kortemme/krivacic/software/anaconda3/lib/'
        process = Popen([rifgen, '@', flags], stdout=PIPE, stderr=PIPE,
                env=myenv)
        stdout, stderr = process.communicate()
        exit_code = process.wait()
        if exit_code != 0:
            print(stdout)
            print(stderr)
        output = subprocess.check_output([rifgen, '@', flags])
        print(output)


if __name__=='__main__':
    main()
