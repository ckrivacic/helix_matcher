'''
Script for submitting matching jobs on the cluster.
Usage:
    submit_matcher.py <pdb> [options]
options:
    --max-runtime=TIME  Maximum amount of time the cluster will allow
    the job to run  [default:12:00:00]
    --max-memory=MEM  Maximum memory allocated to the job  [default: 6G]
    --tasks=#, -j  Number of tasks to split the job into  [default: 500]
    --out, -o  Where to save results  [default: results/]
'''
from klab import cluster, process
import docopt
import sys, os, importlib, shutil

def main():
    
    args = docopt.docopt(__doc__)

    max_runtime = args['--max-runtime']
    max_memory = args['--max-memory']
    ntask = args['--tasks']

    qsub_command = 'qsub', '-h', '-cwd'
    qsub_command += '-o', args['--log']
    qsub_command += '-e', args['--log']
    qsub_command += '-t', '1-{0}'.format(ntask)
    qsub_command += '-l', 'h_rt={}'.format(max_runtime)
    qsub_command += '-l', 'mem_free={}'.format(max_memory)
    qsub_command += '-b', 'y'
    qsub_command += 'N', 'helix_matcher'
    qsub_command += os.environ['ROSEASY_PYTHON'],
    qsub_command += 'matcher.py',
    qsub_command += args['<pdb>'],
    qsub_command += '-j', args['--tasks']
    qsub_command += '-o', args['--out']

    status = process.check_output(qsub_command).decode('utf-8')
    status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
    status_match = status_pattern.match(status)

    if not status_match:
        print(status)
        sys.exit()

    job_id = status_match.group(1)

    qrls_command = 'qrls', job_id
    process.check_output(qrls_command)
    print(status, end=' ')

if __name__=='__main__':
    main()
