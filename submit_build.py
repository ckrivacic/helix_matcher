import glob
from klab import process

def find_nstruct(pdb_prefix):
    file_list = glob.glob(pdb_prefix + '*/*.ent.gz')
    return len(file_list)


pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
nstruct = int(find_nstruct(pdb_prefix) / 100)
max_runtime = '10:00:00'
max_memory = '6G'

qsub_command = 'qsub', '-h', '-cwd'
qsub_command += '-o', 'logs/'
qsub_command += '-e', 'logs/'
qsub_command += '-t', '1-{}'.format(nstruct)
qsub_command += '-l', 'h_rt={}'.format(max_runtime)
qsub_command += '-l', 'mem_free={}'.format(max_memory)
qsub_command += '-b', 'y'
qsub_command += '-N', 'helix_scan'
qsub_command += '/wynton/home/kortemme/krivacic/software/anaconda3/bin/python',
qsub_command += 'build_database.py',

status = process.check_output(qsub_command).decode('utf-8')
status_pattern = re.compile(r'Your job-array (\d+).[0-9:-]+ \(".*"\) has been submitted')
status_match = status_pattern.match(status)

if not status_match:
    print(status)
    sys.exit()

qrls_command = 'qrls', job_id
process.check_output(qrls_command)
print(status, end=' ')
