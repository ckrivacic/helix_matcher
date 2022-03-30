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
'''

import docopt
from helix.matching.scan_helices import PoseScanner
import helix.workspace as ws
import gzip
import pandas as pd
import os, sys
from pathlib import Path
from pyrosetta import init
from pyrosetta import pose_from_file
from datetime import date


def main():
    args = docopt.docopt(__doc__)
    print(args)
    df = pd.DataFrame()
    init('-ignore_unrecognized_res')

    if 'SGE_TASK_ID' in os.environ:
        idx = int(os.environ['SGE_TASK_ID']) - 1
    else:
        idx = int(args['--task']) - 1

    ntasks = int(args['--ntasks']) 
    start = idx * ntasks
    stop = idx * ntasks + ntasks - 1
    print('START: {}'.format(start))
    print('STOP: {}'.format(stop))
    folder = args['<pdbfolder>']
    workspace = ws.workspace_from_dir(args['<workspace>'])

    exts = ('.pdb', '.pdb.gz')
    if args['--recursive']:
        files = [str(p) for p in Path(folder).rglob('*') if
                p.name.endswith(exts)]
    else:
        files = [str(p) for p in Path(folder).iterdir() if
                p.name.endswith(exts)]

    print(files)

    for f in files:
        pose = pose_from_file(f)
        scanner = PoseScanner(pose)
        posename = os.path.basename(f).split('.')[0]
        path = os.path.relpath(f, start=workspace.root_dir)
        helices = pd.DataFrame(
                scanner.scan_pose_helices(name=posename,
                    split_chains=args['--split-chains'], path=path)
                )
        df = pd.concat([df, helices], ignore_index=True)

    today = date.today()
    outdir = os.path.join(workspace.project_params_dir, 'database')
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir,
            'helixdf_custom_{0}.'.format(today.strftime("%m-%d-%Y")))
    df.to_pickle(out + 'pkl')
    # Comment this out later - no need for csv in final
    # df.to_csv(out + 'csv')


if __name__=='__main__':
    main()
