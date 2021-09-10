"""
Run PatchMAN to generate potential interacting scaffolds.

Usage:
    helix patchman <workspace> [options]

Options:
    --local, -l  Run locally
    --target=PDB, -t  Only run for a sepcific target if multiple exist
    --make-dirs  Just make the directories and stop.
    --test-run  Mark as test run. Does nothing for now.
    --clear, -o  Overwrite a prevoius run. Gets rid of docked outputs,
    log files, and job info files.
"""
import helix.workspace as ws
import os
import docopt


def main():
    args = docopt.docopt(__doc__)
    sizes = [10, 14, 21, 28]
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..', 'patchman', 'run_patchman.py'
            )
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targets = args['--target']
    else:
        targets = workspace.targets

    for target in targets:
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        inputs = rif_workspace.unclaimed_inputs
