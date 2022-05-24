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
    --delete  Delete non-designed structures
    --buns-penalty  Include a penalty for buried unsat hbonds
    --prune-buns  Prune rotamers that cause BUNS
    --max-memory=GB  How much memory to allocate  [default: 6G]
    --max-runtime=HH:MM:SS  How long to allocate the CPUs  [default: 36:00:00]
    --designs-per-task=INT  How many designs per task  [default: 60]
    --keep-good-rotamers  Keep rotamers on docked helix that beat the
        average score for its environment in the PDB as long as it does not
        cause buried unsatisfied hbonds.
    --align-thresh=PERCENT  When the match has above this threshold
        sequence identity to the patch, add a favornative bonus.
        [default: 101]
    --special-rot  If using keep-good-rotamers, treat them as special
        rotamers with a score bonus instead of fixing them.
    --special-rot-weight=FLOAT  How much to weigh special rotamers
        [default: -1.5]
    --benchmark  Run the design benchmark (figure out what options to
        pass based on folder name)
    --suffix=STR  Add a suffix to saved designs
    --nocst  Don't constrain during fastdesign
    --ramp-cst  Ramp down constraints
    --upweight-interface=WEIGHT  Upweight interface by this weight
    --hold  Submit the jobs with a hold
    --taskrange=TASKS  Run a range of task #s (local only)
    --subprocess  Run the (local) jobs in a background process
    --special-res  As in special_rot, but with a residue type constraint instead of rotamer.
"""
import helix.workspace as ws
from helix import big_jobs
import os
import math
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
            '..', 'patchman', '03_design_patchman.py'
            )
    if not os.path.exists(script_path):
        raise Exception("Error: {} does not exist.".format(script_path))
    if args['--target']:
        targetlist = args['--target'].split(',')
        targets = [workspace.target_rifdock_path(x) for x in targetlist]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        if not os.path.exists(rif_workspace.log_dir):
            os.makedirs(rif_workspace.log_dir, exist_ok=True)

        if args['--clear']:
            rif_workspace.clear_patchman_designs()

        initial_inputs = sorted(glob.glob(
            os.path.join(rif_workspace.focus_dir, 'patch_*',
                workspace.scaffold_prefix + '*', 'docked_full', '*.pdb.gz')
            ))
        if args['--suffix']:
            inputs = []
            for pdb in initial_inputs:
                if not args['--suffix'] in os.path.basename(pdb):
                    inputs.append(pdb)
        else:
            inputs = initial_inputs

        # ntasks = len(inputs)
        des_per_task = int(args['--designs-per-task'])
        ntasks = math.ceil(len(inputs) / des_per_task)

        cmd = workspace.python_path, script_path
        cmd += rif_workspace.focus_dir,
        # cmd += '--align-thresh', args['--align-thresh']

        if args['--delete']:
            cmd += '--delete',

        if args['--nocst']:
            cmd += '--nocst',

        if args['--buns-penalty']:
            cmd += '--buns-penalty',

        if args['--prune-buns']:
            cmd += '--prune-buns',

        if args['--keep-good-rotamers']:
            cmd += '--keep-good-rotamers',

        if args['--special-rot']:
            cmd += '--special-rot',
        if args['--special-rot-weight'] and not args['--benchmark']:
            cmd += '--special-rot-weight', args['--special-rot-weight']

        if args['--upweight-interface']:
            cmd += '--upweight-interface', args['--upweight-interface']

        if args['--special-res']:
            cmd += '--special-res',

        if args['--task']:
            cmd += '--task', args['--task']
            ntasks = 1

        if args['--taskrange']:
            tasks = []
            for t in range(int(args['--taskrange'].split('-')[0]), int(args['--taskrange'].split('-')[1])):
                tasks.append(t)
            ntasks = len(tasks)

        if args['--suffix']:
            cmd += '--suffix', args['--suffix']

        if args['--ramp-cst']:
            cmd += '--ramp-cst',

        cmd += '--designs-per-task', str(des_per_task)
        
        if args['--benchmark']:
            print('DESIGN TYPE')
            designtype = '_'.join(os.path.basename(target).split('.')[0].split('_')[2:])
            print(designtype)
            if designtype == 'buns_penalty':
                cmd += '--buns-penalty',
            elif designtype == 'buns_penalty_pruned':
                cmd += '--buns-penalty', 
                cmd += '--prune-buns',

            elif designtype == 'residue_lock':
                cmd += '--keep-good-rotamers',

            elif designtype == 'specialres':
                cmd += '--keep-good-rotamers',
                cmd += '--special-res',
            
            elif designtype == 'combined':
                cmd += '--keep-good-rotamers',
                cmd += '--special-res',
                cmd += '--buns-penalty',

            elif designtype == 'combined_nomin':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--buns-penalty',
                cmd += '--nocst',

            elif designtype == 'deleteme':
                cmd += '--align-thresh', '70'

            elif designtype == 'specialrot':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',

            elif designtype == 'specialrot_combined':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--buns-penalty',

            elif designtype == 'special_combined_ramp':
                cmd += '--ramp-cst',
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--buns-penalty',

            elif designtype == 'residue_lock_combined':
                cmd += '--keep-good-rotamers',
                cmd += '--buns-penalty',

            elif designtype == 'residue_lock_combined_ramp':
                cmd += '--keep-good-rotamers',
                cmd += '--buns-penalty',
                cmd += '--ramp-cst',

            elif designtype == 'specialrot_3':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--special-rot-weight', '-3.0'

            elif designtype == 'specialrot_3_combined':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--special-rot-weight', '-3.0'
                cmd += '--buns-penalty',

            elif designtype == 'specialrot_3_combined_ramp':
                cmd += '--keep-good-rotamers',
                cmd += '--special-rot',
                cmd += '--special-rot-weight', '-3.0'
                cmd += '--buns-penalty',
                cmd += '--ramp-cst',


        if args['--local']:
            print('Runinng locally')
            if args['--taskrange']:
                for n in tasks:
                    local_cmd = deepcopy(cmd)
                    local_cmd += '--task', str(n)
                    utils.run_command(local_cmd, background=args['--subprocess'])
            else:
                for n in range(1, ntasks + 1):
                    local_cmd = deepcopy(cmd)
                    if not args['--task']:
                        local_cmd += '--task', str(n)
                    utils.run_command(local_cmd, background=args['--subprocess'])

        else:
            # print('Submitting jobs for {}'.format(target))
            # submit.submit(rif_workspace, cmd, distributor='sge',
                    # make_dirs=args['--make-dirs'],
                    # test_run=args['--test-run'], clear=args['--clear'],)
            script_name='design_patchman'
            print('Submitting jobs for {}'.format(target))
            # submit.submit(rif_workspace, cmd, distributor='sge',
                    # make_dirs=args['--make-dirs'],
                    # test_run=args['--test-run'], clear=args['--clear'],
                    # ntasks=ntasks,
                    # )
            # Call big_jobs.submit directly, so that it doesn't care
            # about unclaimed inputs
            print(args['--max-runtime'])
            print(args['--max-memory'])
            try:
                big_jobs.submit(
                        rif_workspace, cmd,
                        nstruct=ntasks,
                        inputs=inputs,
                        max_runtime=args['--max-runtime'],
                        max_memory=args['--max-memory'],
                        test_run=False,
                        job_name=script_name,
                        create_job_info=True,
                        hold=args['--hold'],
                        )
            except:
                print('{} DID NOT SUBMIT'.format(target))