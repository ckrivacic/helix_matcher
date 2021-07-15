"""
Usage: 
    helix rifdock <workspace> [options]

Options:
    --sge  Running on the cluster?
    --task=NUM  For test runs, just run this task.
    --target=PDB, -t  Only run for a specific target
    --make-dirs  Just make the directories and stop. (Should not need
    this option if you ran prep_rifgen.)
    --test-run  Mark as a test run. For this script this does nothing
    for now.

Workspace should be the root workspace.
"""
from helix import submit
import helix.workspace as ws
from helix import submit
import os
import docopt

def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..', 'matching',
            'matcher.py')
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_targets

    for target in targets:

        print('TARGETS')
        print(targets)
        print('TARGET')
        print(target)
        print('TARGET RIFDOCK PATH')
        print(workspace.target_rifdock_path(target))
        rif_workspace = ws.workspace_from_dir(workspace.target_rifdock_path(target))
        print(type(rif_workspace))

        cmd = workspace.python_path, script_path
        cmd += target,

        if args['--task']:
            cmd += '--task', args['--task']

        print('Submitting jobs for {}'.format(target))
        submit.submit(rif_workspace, cmd, distributor='sge',
                make_dirs=args['--make-dirs'], test_run=args['--test-run'])
