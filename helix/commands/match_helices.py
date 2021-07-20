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
    helix match_helices <workspace> [options]

Options:
    --local, -l  Run locally
    --tasks=NUM, -j  How many tasks to split the run into?
    --verbose, -v  Verbose output
    --database=PATH, -d  Database of relative helix orientations. If not
    provided, will search for database in project_params, then
    standard_params.
    --length, -e  Bin by length (requires a length-binned database)
    --angstroms=NUM, -a  How fine to make the distance bins (requires a
    database with the same binning options). Default set by
    settings.yml.
    --degrees=NUM, -g  How fine to make the angle bins. Defaults set by
    settings.yml.
"""
from helix import submit
import helix.workspace as ws
import docopt
import re
import os
import yaml


# def is_default_settings(workspace):
    # '''Checks if settings are default by looking at which settings.yml
    # file got loaded. Not used anymore, can probably delete.'''
    # return os.path.abspath(os.path.dirname(workspace.settings)) == \
            # os.path.abspath(workspace.standard_params_dir)


def is_default_database(workspace):
    '''Checks if database is default by looking at which folder it's
    in'''
    return os.path.abspath(os.path.dirname(workspace.database_path)) == \
            os.path.abspath(workspace.standard_params_dir)


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

    with open(workspace.settings, 'r') as stream:
        try:
            settings = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    check_overwrite = ['--database', '--length', '--angstroms',
            '--degrees']
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
        if settings['match']['--database'] != 'database/':
            # Database is not default, therefore do not default to
            # project-params.
            database = settings['match']['--database']
            db_origin = 'custom settings'
        elif settings['match']['--database'] == 'database/':
            if not is_default_database(workspace):
                database = workspace.database_path
                db_origin = 'project_params/database'
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


    for target in workspace.targets:
        cmd = workspace.python_path, script_path
        cmd += workspace.target_clusters(target),
        for setting in settings['match']:
            if setting != '--database':
                cmd += setting, settings['match'][setting]
            else:
                cmd += setting, database
        print(cmd)
