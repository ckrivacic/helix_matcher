'''
Usage:
    helix 06_score_matches <workspace> [options]

options:
    --local, -l  Run locally

    --clear, -o  Clear scored outputs (for ALL targets)

    --dataframe=PKL, -d  Path to dataframe of helix vectors  [default: nr_dataframes/final.pkl]

    --plot-alphashape, -p  Show a plot that displays the alphashape of the target protein

    --task=INT, -t  Specify a task number if running locally

    --threshold=NUM  Score over which results will not be saved
    [default: 50]

    --ntasks=INT, -n  How many tasks to split each result dataframe into  [default: 1]

    --target=TAR  Only run for a specific target

    --max-memory=6G  How much memory to allocate to the job?  [default: 4G]

    --max-runtime=10:00:00  How much time to allocate to the job?
    [default: 10:00:00]
'''
from helix import workspace as ws
from helix import big_jobs
from helix.utils import utils
import docopt
import subprocess
import os


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
            '..', 'analysis', 'score_matches.py')
    if not os.path.exists(script_path):
        raise("Error: {} does not exist.".format(script_path))

    if args['--target']:
        targets = [os.path.abspath(args['--target'])]
    else:
        targets = workspace.all_match_workspaces
    for target in targets:
        match_workspace = \
                ws.workspace_from_dir(workspace.target_match_path(target))
        dataframes = match_workspace.outputs
    # else:
    #     dataframes = workspace.all_match_outputs

        if args['--clear']:
            workspace.clear_scores()

        ntasks = int(args['--ntasks']) * len(dataframes)

        cmd = workspace.python_path, script_path
        cmd += match_workspace.focus_dir,
        # if args['--target']:
        #     cmd += '--target', args['--target']
        if args['--plot-alphashape']:
            cmd += '--plot-alphashape',
        if args['--ntasks']:
            cmd += '--ntasks', args['--ntasks']
        if args['--threshold']:
            cmd += '--threshold', args['--threshold']

        if args['--local']:
            if args['--task']:
                cmd += '--task', args['--task']
                utils.run_command(cmd)
            else:
                for i in range(1, ntasks + 1):
                    run_cmd = cmd + ('--task', str(i))
                    utils.run_command(run_cmd)

        else:
            print('Submitting command {}'.format(cmd))
            script_name = 'score_matches'
            big_jobs.submit(match_workspace, cmd,
                    nstruct=ntasks,
                    max_memory=args['--max-memory'],
                    max_runtime=args['--max-runtime'],
                    test_run=False,
                    job_name=script_name,
                    create_job_info=False,)
 

if __name__=='__main__':
    main()
