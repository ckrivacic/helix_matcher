"""
Run FastDesign + FlexPepDock on PatchMAN outputs.

Usage:
    helix design_patchman <workspace> [options]

Options:
    --local, -l  Run locally
    --target=PDB, -t  Only run for a sepcific target if multiple exist
    --make-dirs  Just make the directories and stop.
    --test-run  Mark as test run. Does nothing for now.
    --task=INT  Only run a specific task
    --clear, -o  Overwrite a prevoius run. Gets rid of docked outputs,
    log files, and job info files.
    --flexpepdock  Run flexpepdock at the end of the PatchMAN run
    --relax  Run relax on target
    --delete  Delete non-designed structures
"""
import helix.workspace as ws
import os
import docopt
from helix.utils import utils
from helix import submit
import glob
from copy import deepcopy


def main():
    args = docopt.docopt(__doc__)
    # sizes = [10, 14, 21, 28]
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..', 'patchman', 'design_patchman.py'
            )
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        # inputs = rif_workspace.unclaimed_inputs
        inputs = sorted(glob.glob(
            os.path.join(rif_workspace.focus_dir, 'patch_*',
                workspace.scaffold_prefix + '*', 'docked_full', '*.pdb.gz')
            ))
        ntasks = len(inputs)

        cmd = workspace.python_path, script_path
        cmd += target,

        if args['--flexpepdock']:
            cmd += '--flexpepdock'

        if args['--relax']:
            cmd += '--relax',

        if args['--delete']:
            cmd += '--delete',

        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1

        if args['--local']:
            print('Runinng locally')
            for n in range(1, ntasks + 1):
                local_cmd = deepcopy(cmd)
                if not args['--task']:
                    local_cmd += '--task', str(n)
                utils.run_command(local_cmd)

        else:
            print('Submitting jobs for {}'.format(target))
            submit.submit(rif_workspace, cmd, distributor='sge',
                    make_dirs=args['--make-dirs'],
                    test_run=args['--test-run'], clear=args['--clear'],)
