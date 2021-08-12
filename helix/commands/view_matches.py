'''
Make a PyMOL session for each result, one at a time, in ascending score
order.

Usage:
    helix view_matches <workspace> [options]

Options:
    --target, -t  Only view matches for the given target
'''
from helix.analysis import clash
from helix import workspace as ws
import docopt
import prody
import pandas as pd
import pickle5 as pickle
import pymol
from helix.matching import scan_helices
import networkx as nx
from pyrosetta import pose_from_file
from pyrosetta import init
from helix.utils import numeric
import subprocess
import os

def session_from_graph(results_row, query_df, db_df):

    init()
    def chain_from_name(string):
        chainno = int(string.split('_')[-1])
        chain = chr(ord('@')+chainno)
        return chain

    subgraph = results_row['subgraph']

    query_selstr = ""
    db_selstr = "db and "

    # query_path = os.path.dirname(query_df.loc[0]['path'])

    db_sels = []
    df_row = db_df.loc[subgraph[0][0]]
    if 'path' in df_row:
        pdb = df_row['path']
        pdb_name = os.path.basename(pdb)
        # pymol.cmd.load(pdb, pdb_name)
        dfpose = pose_from_file(pdb)
    else:
        pdb = df_row['name'].split('_')[0]
        pdb_name = pdb
        # pymol.cmd.fetch(pdb, pdb_name)
        dfpose = pose_from_file(pdb + '.cif')
    query_vectors = []
    df_vectors = []
    qobjs = '('

    query_paths = []
    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]
        query_paths.append(query_row['path'])

        try:
            start = dfpose.pdb_info().pose2pdb(df_row['start']).split(' ')[0]
            stop = dfpose.pdb_info().pose2pdb(df_row['stop']).split(' ')[0]
        except:
            start = df_row['start']
            stop = df_row['stop']

        query_vectors.extend(query_row['vector'])
        df_vectors.extend(df_row['vector'])

        query_name = os.path.basename(query_row['path'])[:-7]
        qobjs += query_name + ' or '

        query_selstr += "({} and (resi {}-{})) or ".format(
                query_name,
                query_row['start'], query_row['stop']
                )
        db_selstr = "(resi {}-{} and chain {})".format(
                start, stop,
                df_row['chain']
                )
        db_sels.append(db_selstr)

    transform = numeric.Transformation(df_vectors, query_vectors)
    pymol_transform = get_pymol_transform(transform)

    db_selstr = '('+ db_sels[0]
    for i in range(1, len(db_sels)):
        db_selstr += ' or '
        db_selstr += db_sels[i]
    db_selstr += ')'
    print(db_selstr)

    query_selstr = query_selstr[:-4]
    query_selstr += ' and chain A'

    qobjs = qobjs[:-3] + ')'

    print(query_selstr)
    print(db_selstr)

    script_path = os.path.dirname(os.path.realpath(__file__))
    script_path = os.path.join(script_path, '..', 'analysis',
            'launch_pymol.sho')
    cmd = [script_path]
    cmd.append(';'.join(query_paths))
    cmd.append(pdb)
    cmd.append(query_selstr)
    cmd.append(db_selstr)
    cmd.append(df_row['chain'])
    cmd.extend(pymol_transform)
    cmd.append(qobjs)

    print(cmd)

    subprocess.call(cmd)


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])

    try:
        df = pd.read_pickle(workspace.dataframe_path)
        helices = pd.read_pickle(workspace.all_scaffold_dataframe)
    except:
        with open(workspace.dataframe_path, 'rb') as f:
            df = pickle.load(f)
        with open(workspace.all_scaffold_dataframe, 'rb') as f:
            helices = pickle.load(f)

    for target in workspace.targets:
        scored_outputs = workspace.match_outputs(target)
        results = pd.DataFrame()
        for outfile in scored_outputs:
            try:
                scores = pd.read_pickle(outfile)
            except:
                with open(outfile, 'rb') as f:
                    scores = pickle.load(f)
            results = results.concat([results, scores],
                    ignore_index=True)
        results = results.sort_values(by='total_match_score',
                ascending=True)
        for i in range(0, 100):
            testrow = results.iloc[i]
            session_from_graph(testrow, helices, df, alpha)
