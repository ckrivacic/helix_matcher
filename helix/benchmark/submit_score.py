'''
Usage: submit_score.py [options]

Options:
    --min  Run minimized
    --buried  Get buried residues instead of interface
    --cst-bb  Constrain backbone coordinates
    --cst-sc  Constrain sidechain coordinates
'''
import glob, re, os
from klab import process
from helix.utils import utils
import docopt
import sys

args = docopt.docopt(__doc__)
nstruct = len(utils.get_pdb_list(pdbid=False))
minimize = args['--min']
buried = args['--buried']
cst_bb = args['--cst-bb']
cst_sc = args['--cst-sc']

if minimize:
    num=50
else:
    num=250
ntasks = (nstruct // num) + 1
max_runtime = '24:00:00'
max_memory = '6G'
logdir = 'residue_scores_min/logs/'
os.makedirs(logdir, exist_ok=True)

qsub_command = 'qsub', '-h', '-cwd'
qsub_command += '-o', logdir 
qsub_command += '-e', logdir
qsub_command += '-t', '1-{}'.format(ntasks)
qsub_command += '-l', 'h_rt={}'.format(max_runtime)
qsub_command += '-l', 'mem_free={}'.format(max_memory)
qsub_command += '-b', 'y'
qsub_command += '-N', 'pdb_score'
qsub_command += '/wynton/home/kortemme/krivacic/software/anaconda3/bin/python',
qsub_command += 'score_pdb.py', str(num)
if minimize:
    qsub_command += '--min',
if buried:
    qsub_command += '--buried',
if cst_bb:
    qsub_command += '--cst-bb',
if cst_sc:
    qsub_command += '--cst-sc',

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
