"""
Runs PatchMAN.

Usage:
    patchman.py <workspace> [options]

Options:
    --target=PDB, -t  Only runs on a specific target
"""
import docopt
from helix import workspace as ws


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
