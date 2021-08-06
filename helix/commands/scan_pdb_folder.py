
'''
Create a helix dataframe from a folder containing PDBs.

Usage:
    helix scan_pdb_folder <workspace> <pdbfolder> [options]

Options:
    --recursive, -r  Scan subfolders as well
    --ntasks=NUM, -n  How many tasks to split the job into  [default: 1]
    --task=NUM, -t  Which task for this job (if not using SGE)?
    1-indexed  [default: 1]
    --split-chains, -s  Split the pose into separate chains?  [default: True]
    --local, -l  Run locally (recommended)  [default: True]
    --sge  Run on SGE (for very large folders)
    --clear  Deletes the existing custom database
    --max-memory=6G  How much memory to allocate to the job?  [default: 4G]

    --max-runtime=10:00:00  How much time to allocate to the job?
    [default: 10:00:00]
'''
import docopt
from helix import workspace as ws
from copy import deepcopy
from helix.utils import utils
import os
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
        workspace.clear_database()

    cmd = workspace.python_path, script_path
    cmd += workspace.root_dir, args['<pdbfolder>']
    argpass = ['--ntasks']
    for arg in argpass:
        cmd += arg, str(args[arg])
    if args['--recursive']:
        cmd += '--recursive',
    if args['--split-chains']:
        cmd += '--split-chains'
    ntasks = int(args['--ntasks'])

    if not os.path.exists(workspace.log_dir):
        os.makedirs(workspace.log_dir, exist_ok=True)

    if args['--local']:
        if args['--task']:
            local_cmd = deepcopy(cmd)
            utils.run_command(local_cmd)
        else:
            for n in range(1, ntasks + 1):
                local_cmd= deepcopy(cmd)
                local_cmd += '--task', str(n)
                utils.run_command(local_cmd)

    else:
        script_name = 'scan_helices'
        print('Submitting the following command to SGE:')
        print(' '.join(cmd))
        big_jobs.submit(
                workspace, cmd, nstruct=ntasks,
                max_runtime=args['--max-runtime'],
                max_memory=args['--max-memory'],
                job_name='scan_helices',
                test_run=False,
                create_job_info=False
                )
