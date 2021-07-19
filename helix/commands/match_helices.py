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

    cmd = workspace.python_path, script_path
    for setting in settings['match']:
        cmd += setting, settings['match'][setting]
    print(cmd)
