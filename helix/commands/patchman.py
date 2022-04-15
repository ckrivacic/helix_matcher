"""
Run PatchMAN to generate potential interacting scaffolds.

Usage:
    helix patchman <workspace> [options]

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
    --keep-existing, -k  Keep existing results; only run for patches
        that do not have results.
    --max-runtime TIME     [default: 12:00:00]
        The runtime limit for each job.
    --max-memory MEM       [default: 2G]
        The memory limit for each job.
    --hold, -h  Keep the job on hold after submitting (to be released
        manually later)
    --taskrange=TASKS  Run a range of task #s (local only)
    --subprocess  Run the (local) jobs in a background process
"""
import helix.workspace as ws
import os
import docopt
from helix.utils import utils
from helix import big_jobs
from helix import submit
from copy import deepcopy


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
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        rif_workspace.make_dirs()
        if args['--make-dirs']:
            continue
        if args['--clear']:
            rif_workspace.clear_outputs
        if args['--keep-existing']:
            inputs = rif_workspace.unfinished_inputs
            print('Unfinished inputs')
            print(inputs)
        else:
            inputs = rif_workspace.unclaimed_inputs
        ntasks = len(inputs)
        if ntasks == 0:
            print('No inputs for target {}'.format(target))
            continue

        cmd = workspace.python_path, script_path
        cmd += target,

        if args['--flexpepdock']:
            cmd += '--flexpepdock'

        if args['--relax']:
            cmd += '--relax'

        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1

        if args['--taskrange']:
            tasks = []
            for t in range(int(args['--taskrange'].split('-')[0]), int(args['--taskrange'].split('-')[1])):
                tasks.append(t)
            ntasks = len(tasks)

        if args['--local']:
            print('Runinng locally')
            if args['--taskrange']:
                for n in tasks:
                    local_cmd = deepcopy(cmd)
                    local_cmd += '--task', str(n)
                    utils.run_command(local_cmd, background=args['--subprocess'], logdir=rif_workspace.log_dir,
                                      log_prefix='task_{}'.format(n))
            else:
                for n in range(1, ntasks + 1):
                    local_cmd = deepcopy(cmd)
                    if not args['--task']:
                        local_cmd += '--task', str(n)
                    utils.run_command(local_cmd, background=args['--subprocess'])

        else:
            print('Submitting jobs for {}'.format(target))
            big_jobs.submit(
                    rif_workspace, cmd, nstruct=ntasks, 
                    max_runtime=args['--max-runtime'],
                    max_memory=args['--max-memory'],
                    job_name='PatchMAN',
                    inputs=inputs,
                    hold=args['--hold'],
                    )
            # submit.submit(rif_workspace, cmd, distributor='sge',
                    # make_dirs=args['--make-dirs'],
                    # test_run=args['--test-run'], clear=args['--clear'],
                    # max_runtime='24:00:00', max_memory='6G',
                    # inputs=inputs)
