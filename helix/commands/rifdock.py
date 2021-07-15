"""
Usage: 
    helix rifdock <workspace> [options]

Options:
    --sge  Running on the cluster?
    --task=NUM  For test runs, just run this task.
    --target=PDB, -t  Only run for a specific target

Workspace should be the root workspace.
"""
from helix import submit
import helix.workspace as ws
from helix import submit
import docopt

def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.realpath(__file__),
            '..', 'matching',
            'matcher.py')
    if not os.path.exists(script_path):
        raise("Error: matcher.py does not exist.")
    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_targets

    if args['--database']:
        db = args['--database']
    else:
        db = workspace.database_path

    for target in targets:

        rif_workspace = ws.workspace_from_dir(workspace.target_rifdock_path(target))

        cmd = workspace.python_path, script_path
        cmd += target
        cmd += '--database', db 

        if args['--task']:
            cmd += '--task', args['--task']

        print('Submitting jobs for {}'.format(target))
        submit.submit(rif_workspace, cmd, distributor='sge')
