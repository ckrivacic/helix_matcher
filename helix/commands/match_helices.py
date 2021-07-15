"""
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
from helix import submit
import docopt

def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    script_path = os.path.join(
            os.path.realpath(__file__),
            '..', 'matching',
            'matcher.py')
    if not os.path.exists(script_path):
        raise("Error: matcher.py does not exist.")

    if args['--database']:
        db = args['--database']
    else:
        db = workspace.database_path

    cmd = workspace.python_path, script_path
    cmd += '--database', db 
