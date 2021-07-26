
'''
Create a helix dataframe from a folder containing PDBs.

Usage:
    scan_pdb_folder.py <workspace> <pdbfolder> [options]

Options:
    --recursive, -r  Scan subfolders as well
    --ntasks=NUM, -n  How many tasks to split the job into  [default: 1]
    --task=NUM, -t  Which task for this job (if not using SGE)?
    1-indexed  [default: 1]
    --split-chains, -s  Split the pose into separate chains?  [default: True]
    --local, -l  Run locally (recommended)  [default: True]
    --sge  Run on SGE (for very large folders)
'''
import docopt
from helix import workspace as ws
from copy import deepcopy
from helix.utils import utils
from helix import big_jobs


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    pdffolder = args['<pdbfolder>']

    script_path = \
            os.path.join(os.path.abspath(
                os.path.dirname(os.path.realpath(__file__))), '..',
            'matching', 'scan_pdb_folder.py')
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))

    if args['--sge']:
        args['--local'] = False
    if args['--clear']:
        workspace.clear_database

    cmd = workspace.python_path, script_path
    argpass = ['--recursive', '--ntasks', '--split-chains']


    if args['--local']:
        if args['--task']:
            local_cmd = deepcopy(cmd)
            for arg in argpass:
                local_cmd += arg, args[arg]
            utils.run_command(local_cmd)
        else:
            for n in range(1, args['--ntasks'] + 1):
                local_cmd= deepcopy(cmd)
                local_cmd += '--task', str(n)
                utils.run_command(local_cmd)

    else:
        script_name = 'scan_helices'
        print('Submitting the following command to SGE:')
        print(' '.join(cmd))
        big_jobs.submit(
                workspace, cmd, nustruct=ntasks,
                max_runtime=args['--max-runtime'],
                max_memory=args['--max-memory'],
                test_run=False,
                create_job_info=False
                )
