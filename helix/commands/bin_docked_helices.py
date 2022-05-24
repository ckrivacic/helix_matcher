"""
Calculates helical vectors and bins them according to the Angstrom and degree bins provided in the workspace's settings
or via the options given to this command.

Usage: 
    helix 04_match <workspace> [options]

Options:
    --local, -l  Run locally

    --verbose, -v  Verbose output

    --clear, -o  Clear outputs (logs + dataframes) before starting run

    --make-dirs  Make the necessary directories and then stop

    --database=PATH, -d  Database of relative helix orientations. If not
    provided, will search for database in project_params, then
    standard_params.

    --ntasks=NUM, -n  How many sub-tasks to break up each matching run
    into. The total number of tasks will be (# targets)(# dataframes in
    database)(tasks), but each job will have (# dataframes)*(tasks)
    [default: 1]

    --length, -e  Bin by length (requires a length-binned database)

    --angstroms=NUM, -a  How fine to make the distance bins (requires a
    database with the same binning options). Default set by
    settings.yml.

    --degrees=NUM, -g  How fine to make the angle bins. Defaults set by
    settings.yml.

    --target=PDBPATH, -t  Only run matcher on the given target.

    --scaffold=SCAFFOLD  Only run matcher for a given docked scaffold
    (meaning helix length for now).

    --max-memory=6G  How much memory to allocate to the job?  [default: 4G]

    --max-runtime=10:00:00  How much time to allocate to the job?
    [default: 10:00:00]
"""
from helix import submit
from helix.utils import utils
from helix import big_jobs
import helix.workspace as ws
import docopt
import re
import os
import yaml


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            '..', 'matching',
            'matcher.py')
    if not os.path.exists(script_path):
        print(script_path)
        raise("Error: matcher.py does not exist.")

    if args['--database']:
        db = args['--database']
    else:
        db = workspace.database_path

    # with open(workspace.settings, 'r') as stream:
        # try:
            # settings = yaml.safe_load(stream)
        # except yaml.YAMLError as exc:
            # print(exc)

    settings = workspace.settings
    check_overwrite = ['--length', '--angstroms',
            '--degrees']
    # For these settings, use command-line arguments over inputs from
    # settings file.
    for setting in check_overwrite:
        if args[setting]:
            settings['match'][setting] = args[setting]

    db_origin = 'default'
    if not args['--database']:
        # This is the default database path. If it has not been
        # modified, make sure to append the subdirectory which
        # corresponds to the actual database. This code may need to be
        # modified in the future depending on what other database types
        # end up in the default package.
        if not workspace.is_default_database(settings['match']['--database']):
            # Database is not default, therefore do not default to
            # project-params.
            database = settings['match']['--database']
            db_origin = 'custom settings'
        else:
            # If the settings database is default, see if there is a
            # non-default database in the workspace (usually
            # project_params)
            if not workspace.is_default_database(workspace.database_path):
                database = workspace.database_path
                db_origin = 'project_params/database'
            # Otherwise, use the default settings database.
            else:
                database = settings['match']['--database']

    else:
        database = args['--database']
        db_origin = 'command argument'

    dbexp = 'bins_.*A_.*D'
    match = False
    for subfolder in os.listdir(database):
        if re.match(dbexp, subfolder):
            match = True
            break

    if not match:
        if args['--length']:
            db_subdir = 'length'
        else:
            db_subdir = 'standard'
        database = os.path.join(database, db_subdir)

    if not os.path.exists(database):
        sys.exit("Could not find database at {}. Make sure your database "\
                "path is correct. Database determined via "\
                "{}.".format(database, db_origin))

    else:
        print("Matching using database at {}, determined via "\
                "{}".format(database, db_origin))

    if args['--target']:
        targets = [os.path.abspath(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    for target in targets:
        match_workspace = ws.MatchWorkspace(workspace.root_dir, target)
        match_workspace.make_dirs()

        if args['--clear']:
            match_workspace.clear_outputs()

        if args['--make-dirs']:
            continue

        cmd = match_workspace.python_path, script_path
        cmd += 'bin', match_workspace.scaffold_dataframe
        # cmd += match_workspace.target_clusters(target),
        for setting in settings['match']:
            if setting != '--database':
                cmd += setting, settings['match'][setting]
            else:
                cmd += setting, database
        if args['--ntasks']:
            ntasks = int(args['--ntasks']) * workspace.n_bin_pickles
            cmd += '--tasks', args['--ntasks']

        if args['--scaffold']:
            cmd += '--scaffold', args['--scaffold']

        if args['--local']:
            cmd += '--local',
            utils.run_command(cmd)
            continue
        else:
            script_name = 'matcher'
            print('Submitting jobs for {}'.format(target))
            big_jobs.submit(match_workspace, cmd,
                    nstruct=ntasks,
                    max_memory=args['--max-memory'],
                    max_runtime=args['--max-runtime'],
                    test_run=False,
                    job_name=script_name,
                    create_job_info=False,)
