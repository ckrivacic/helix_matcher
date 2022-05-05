"""
Run the matching algorithm. Requires a database containing the relative
orientations of helices in a variety of scaffold proteins (default
database contains all proteins in the nonredundant PDB). The relative
orientations in this database must be binned with the same intervals as
are given here; if they are not, you will have to start from a
vectorized helices dataframe and bin them according to the settings you
want (use the "helix bin" command NOT IMPLEMENTED).

Settings for matching can be set using the settings.yml file. Default
settings are in standard_params/settings.yml. If you place a
settings.yml file in project_params, those will be loaded instead.
(Note: If you do this, make sure you include values for all settings
present in the default file. The suggested method is to copy the default
file over and edit it from there.)
Additionally, you can supply options to this submission script, and
those will overwrite any settings from the yaml file(s).

Database defaults are prioritized as follows:
    1. A database passed to this command via --database option
    2. A database OTHER THAN the default database path passed via
    settings.yml
    3. A database placed in project_params/
    4. The default database

If the database does not have a bins_.*A_.*D regex pattern in its
subfolders, this script will try to look one level down for
'standard' and 'length' subfolders. If not found, then you should fix
your database path.

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

    --min-dist=FLOAT  Minimum distance for two helices to have their relative orientations saved.
    Useful for query dataframes where you have lots of helices.

    --max-dist=FLOAT  Maximum distance before two helices do not have their relative orientations saved.

    --bin-tasks=INT  Split binning into this many tasks. This will also turn this into a bin-only run, i.e.
    no matching. Run this command again to match using the created bins.
"""
from helix import submit
from helix.utils import utils
from helix import big_jobs
from copy import deepcopy
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
        if args['--bin-tasks']:
            cmd += 'bin_query', match_workspace.focus_dir
            cmd += '--bin-tasks', args['--bin-tasks']
        else:
            cmd += 'match', match_workspace.focus_dir
        # cmd += match_workspace.target_clusters(target),
        for setting in settings['match']:
            if setting != '--database':
                cmd += setting, settings['match'][setting]
            else:
                cmd += setting, database
        num_rel_dataframes = len(match_workspace.relative_orientation_dataframes)
        if args['--ntasks'] and not args['--bin-tasks']:
            ntasks = int(args['--ntasks']) * workspace.n_bin_pickles# * num_rel_dataframes
            cmd += '--tasks', args['--ntasks']
        if args['--bin-tasks']:
            ntasks = int(args['--bin-tasks'])

        if args['--scaffold']:
            cmd += '--scaffold', args['--scaffold']

        if args['--min-dist']:
            cmd += '--min-dist', args['--min-dist']
        if args['--max-dist']:
            cmd += '--max-dist', args['--max-dist']

        if args['--verbose']:
            cmd += '--verbose',

        if args['--local']:
            cmd += '--local',
            if args['--bin-tasks']:
                for i in range(1, int(args['--bin-tasks']) + 1):
                    cmd_copy = deepcopy(cmd)
                    cmd_copy += '--bin-task-id', str(i)
                    utils.run_command(cmd_copy)
            else:
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
