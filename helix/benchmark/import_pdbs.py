"""
Usage:
    import_pdbs.py <workspace> <dataframe>
"""
import pandas as pd
import os
import docopt
import yaml
from helix import workspace as ws


def link_pdb(workspace, pdbid, chain):
    pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
    path = os.path.join(pdb_prefix, pdbid[1:3],
            'pdb{}.ent.gz'.format(pdbid))
    os.symlink(path, os.path.join(workspace.target_dir,
        '{}_{}.pdb.gz'.format(pdbid, chain)))


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    try:
        df = pd.read_pickle(args['<dataframe>'])
    except:
        import pickle5 as pickle
        with open (args['<dataframe>'], 'rb') as f:
            df = pickle.load(f)

    chainmap = {}

    for idx, row in df.iterrows():
        pdbid = row['name']
        chain = row['target']
        link_pdb(workspace, pdbid, chain)
        chainmap['{}_{}'.format(pdbid, chain)] = chain

    chainmap_path = os.path.join(
            workspace.project_params_dir,
            'chainmap.yml'
            )

    stream = file(chainmap_path, 'w')
    yaml.dump(chainmap, stream)
