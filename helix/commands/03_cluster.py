"""
Cluster RIFDock outputs. The lowest scoring docking output is chosen as
the representative of each cluster and copied into the
cluster_representatives directory.

Usage: 
    helix 03_cluster <workspace> [options]

Options:
    --local, -l  Run locally, each target/docked scaffold in sequence.
    Cluster runs should not be necessary.
    --task=NUM  For test runs, just run this task.
    --target=PDB, -t  Only run for a specific target
    --make-dirs  Just make the directories and stop. (Should not need
    this option if you ran prep_rifgen.)
    --test-run  Mark as a test run. For this script this does nothing
    for now.
    --clear, -o  Overwrite a previous run. Gets rid of docked outputs,
    log files, and job info files.
    --max-memory=4G  How much memory to allocate to the job?  [default: 4G]
    --max-runtime=4:00:00  How much time to allocate to the job?
    [default: 4:00:00]

Workspace should be the root workspace.
"""
from helix import big_jobs
import helix.workspace as ws
from helix import submit
from helix.utils import utils
import os, sys
import docopt
from copy import deepcopy

def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..', 'rifdock',
            'cluster_sge.py')
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:

        rif_workspace = ws.workspace_from_dir(target)
        inputs = 1
        scaffolds = workspace.helices
        ntasks = len(scaffolds)

        cmd = workspace.python_path, script_path
        cmd += target,

        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1
            if args['--local']:
                utils.run_command(cmd)
                sys.exit()

        if args['--local']:
            for n in range(1, ntasks + 1):
                local_cmd = deepcopy(cmd)
                if not args['--task']:
                    local_cmd += '--task', str(n),
                utils.run_command(local_cmd)
        else:
            script_name='cluster'
            print('Submitting jobs for {}'.format(target))
            # submit.submit(rif_workspace, cmd, distributor='sge',
                    # make_dirs=args['--make-dirs'],
                    # test_run=args['--test-run'], clear=args['--clear'],
                    # ntasks=ntasks,
                    # )
            if args['--clear']:
                rif_workspace.clear_cluster_outputs()
            print('Submitting the following command to SGE:')
            print(' '.join(cmd))
            # Call big_jobs.submit directly, so that it doesn't care
            # about unclaimed inputs
            big_jobs.submit(
                    rif_workspace, cmd,
                    nstruct=ntasks,
                    max_runtime=args['--max-runtime'],
                    max_memory=args['--max-memory'],
                    test_run=False,
                    job_name=script_name,
                    create_job_info=False,
                    )
