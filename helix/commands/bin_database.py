"""
Bin a database.

Usage:
    helix bin_database <workspace> [options]

Options:
    --database=PATH, -d  If database is not in the workspace, provide a
    path.

    --local, -l  run locally

    --verbose, -v  Verbose output

    --clear  Clear the database before running

    --ntasks=NUM, -n  How many sub-tasks to break the binning process
    into?

    --length, -e  Bin by length

    --angstroms=NUM, -a  How fine to make the distance bins (requires a
    database with the same binning options). Default set by
    settings.yml.

    --degrees=NUM, -g  How fine to make the angle bins. Defaults set by
    settings.yml.

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
import sys
import glob
import yaml


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.abspath(os.path.dirname(os.path.realpath(__file__))),
            '..', 'matching',
            'matcher.py')

    if args['--database']:
        db = args['--database']
    else:
        db = workspace.database_path

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
            print(settings['match']['--database'])
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

    print('Database at {} being used based on {}'.format(database, db_origin))

    if args['--length']:
        out = os.path.join(database, 'length')
    else:
        out = os.path.join(database, 'standard')

    os.makedirs(out, exist_ok=True)
    picklepath = glob.glob(os.path.join(database, 'helixdf*.pkl'))
    if len(picklepath) > 1:
        print("Multiple helix dataframes found in {0}! Please consolodate them "\
                "or remove extraneous pkl files beginning with "\
                "'helixdf'.".format(database))
        sys.exit()
    elif len(picklepath) == 0:
        print("No helix dataframe ('helixdf*.pkl') found in provided "\
                "database path ({0}). Please runs 'helix scan' command "\
                "before binning.".format(database))
        sys.exit()

    helixdf = picklepath[0]

    if args['--clear']:
        workspace.clear_database()

    cmd = workspace.python_path, script_path
    cmd += 'bin', helixdf
    for setting in settings['match']:
        if setting != '--database':
            cmd += setting, settings['match'][setting]
        else:
            cmd += setting, database

    if args['--ntasks']:
        cmd += '--tasks', args['--tasks']

    if args['--local']:
        cmd += '--local',
        utils.run_command(cmd)
    
    else:
        script_name='helixbin'
        big_jobs.submit(workspace, cmd, 
                nstruct=args['--ntasks'],
                max_memory=args['--max-memory'],
                max_runtime=args['--max-runtime'],
                test_run=False,
                job_name=script_name,
                create_job_info=False)
