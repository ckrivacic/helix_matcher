"""
Run RIFGen and RIFDock. Must be run after 01_prep_rifgen. By default,
this attempts to submit these jobs to an SGE job distributor. You can
run it locally using the -l command. For SGE runs, each target is given
its own job, and each patch on the target surface has its own task
within that job. 

Due to the size of the rotamer interaction fields,
RIFGen outputs are stored in a temp directory, which on Wynton is
automatically deleted after job completion, and only the final
RIFDock outputs are transferred back into the workspace.

Usage: 
    helix 02_rifdock <workspace> [options]

Options:
    --local, -l  Run locally, each patch/target in sequence
    --task=NUM  For test runs, just run this task.
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
from helix.utils import utils
import os
import docopt
from copy import deepcopy

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

        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        inputs = rif_workspace.unclaimed_inputs
        ntasks = len(inputs)

        cmd = workspace.python_path, script_path
        cmd += target,

        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1

        if args['--local']:
            local_cmd = deepcopy(cmd)
            for n in range(1, ntasks + 1):
                if not args['--task']:
                    local_cmd += '--task', str(n),
                utils.run_command(local_cmd)
        else:
            cmd += '--sge',
        print('Submitting jobs for {}'.format(target))
        submit.submit(rif_workspace, cmd, distributor='sge',
                make_dirs=args['--make-dirs'],
                test_run=args['--test-run'], clear=args['--clear'],
                )
