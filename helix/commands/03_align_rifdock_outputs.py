"""
Usage: 
    helix 03_align_rifdock_outputs <workspace> [options]

Options:
    --local, -l  Run serially on a local computer instead of in parallel
    on the cluster.
    --target=PDB, -t  Only run for a specific target
    --make-dirs  Just make the directories and stop. (Should not need
    this option if you ran prep_rifgen.)
    --test-run  Mark as a test run. For this script this does nothing
    for now.
    --clear, -o  Overwrite a previous run. Gets rid of docked outputs,
    log files, and job info files.

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
            '..', 'rifdock',
            'run_sge.py')
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:

        rif_workspace = ws.workspace_from_dir(target)

        cmd = workspace.python_path, script_path
        cmd += target,

        ntasks = None
        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1

        print('Submitting jobs for {}'.format(target))
        submit.submit(rif_workspace, cmd, distributor='sge',
                make_dirs=args['--make-dirs'],
                test_run=args['--test-run'], clear=args['--clear'])
