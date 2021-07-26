
'''
Create a helix dataframe from a folder containing PDBs.

Usage:
    scan_pdb_folder.py <workspace> <pdbfolder> [options]

Options:
    --recursive, -r  Scan subfolders as well
    --ntasks=NUM, -n  How many tasks to split the job into  [default: 1]
    --task=NUM, -t  Which task for this job (if not using SGE)?
    1-indexed  [default: 1]
    --split-chains, -s  Split the pose into separate chains?  [default: True]
    --local, -l  Run locally (recommended)  [default: True]
    --sge  Run on SGE (for very large folders)
'''
import docopt
from helix import workspace as ws


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    pdffolder = args['<pdbfolder>']

    if args['--sge']:
        args['--local'] = False
        script_name = 'scan_helices'
        if args['--clear']:
            workspace.clear_database
