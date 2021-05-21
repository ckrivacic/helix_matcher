import pandas as pd
import subprocess
import sys, os, re
import pickle
from klab.rosetta import input_files
from roseasy import pipeline
from roseasy.big_jobs import run_command

'''
Usage:
    python3 submit_cluster.py <folder>
'''


def submit_qsub(parent, folders):
    python_path = os.environ['ROSEASY_PYTHON'] 
    max_runtime = '24:00:00'
    max_memory = '6G'
    # ntasks = ntasks * len(tasks)
    script_path = os.path.abspath(
            'cluster.py'
            )
    ntasks = len(folders)
    logdir = os.path.join(parent, 'cluster', 'logs')
    os.makedirs(logdir, exist_ok=True)
    qsub_command = 'qsub', '-h', '-cwd',
    qsub_command += '-o', logdir,
    qsub_command += '-t', '1-{}'.format(ntasks),
    qsub_command += '-l', 'h_rt={}'.format(max_runtime),
    qsub_command += '-l', 'mem_free={}'.format(max_memory),
    qsub_command += '-b', 'y',
    qsub_command += '-N', 'cluster',
    qsub_command += workspace.python_path,
    qsub_command += script_path,
    qsub_command += parent,

    status = subprocess.check_output(qsub_command).decode('utf-8')
    status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
    status_match = status_pattern.match(status)

    if not status_match:
        print(status)
        sys.exit()

    # Figure out the job id, then make a params file specifically for it.

    job_id = status_match.group(1)

    # with open(workspace.job_info_path(job_id), 'w') as file:
        # json.dump(params, file)

    # Release the hold on the job.

    qrls_command = 'qrls', job_id
    subprocess.check_output(qrls_command)
    print(status, end=' ')


def main():
    import glob
    parent = sys.argv[1]
    folders = glob.glob('{}/*_output')
    submit_qsub(parent, folders)


if __name__=='__main__':
    main()
