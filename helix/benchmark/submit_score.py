import glob, re, os
from klab import process
from helix.utils import utils

nstruct = len(utils.get_pdb_list(pdbid=False))
num=1000
ntasks = nstruct / num + 1
max_runtime = '24:00:00'
max_memory = '6G'
logdir = 'residue_scores/logs/'
os.makedirs(logdir, exist_ok=True)

qsub_command = 'qsub', '-h', '-cwd'
qsub_command += '-o', logdir 
qsub_command += '-e', logdir
qsub_command += '-t', '1-{}'.format(nstruct)
qsub_command += '-l', 'h_rt={}'.format(max_runtime)
qsub_command += '-l', 'mem_free={}'.format(max_memory)
qsub_command += '-b', 'y'
qsub_command += '-N', 'pdb_score'
qsub_command += '/wynton/home/kortemme/krivacic/software/anaconda3/bin/python',
qsub_command += 'score_pdb.py', str(num)

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
