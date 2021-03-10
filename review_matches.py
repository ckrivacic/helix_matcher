'''
Usage:
    review_matches.py <results> [options]

options:
    --dataframe=PKL, -d
        Path to dataframe of helix vectors  [default: nr_dataframes/final.pkl]
'''

import docopt
import pandas as pd
import pymol
import scan_helices
import networkx as nx
from pyrosetta import pose_from_file
from pyrosetta import init
import numeric
import subprocess
import os


def max_subgraph(graph):
    max_subgraph_len = 0
    max_subgraph = None
    subgraphs = []
    for f in nx.find_cliques(graph):
        if len(f) > max_subgraph_len:
            subgraphs = []
            max_subgraph_len = len(f)
            subgraphs.append(f)
        elif len(f) == max_subgraph_len:
            subgraphs.append(f)

    return subgraphs

def session_from_graph(graph, query_df, db_df):

    def chain_from_name(string):
        chainno = int(string.split('_')[-1])
        chain = chr(ord('@')+chainno)
        return chain

    subgraph = max_subgraph(graph)[0]
    print(subgraph)
    query_selstr = ""
    db_selstr = "db and "

    for name, group in query_df.groupby('path'):
        print('Loading {} into pymol'.format(
            name
            ))
        pymol.cmd.load(name, os.path.basename(name[:-7]))

    query_path = os.path.dirname(query_df.loc[0]['path'])

    db_sels = []
    df_row = db_df.loc[subgraph[0][0]]
    pdb = df_row['name'].split('_')[0]
    print('fetching {}'.format(
        df_row['name'].split('_')[0]
        ))
    pymol.cmd.fetch(df_row['name'].split('_')[0], 'db')
    dfpose = pose_from_file(df_row['name'].split('_')[0] + '.cif')
    query_vectors = []
    df_vectors = []
    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]

        # start = dfpose.pdb_info().pose2pdb(df_row['start']).split(' ')[0]
        # stop = dfpose.pdb_info().pose2pdb(df_row['stop']).split(' ')[0]
        start = df_row['start']
        stop = df_row['stop']
        
        query_vectors.extend(query_row['vector'])
        df_vectors.extend(df_row['vector'])

        query_selstr += "({} and (resi {}-{})) or ".format(
                os.path.basename(query_row['path'])[:-7],
                query_row['start'], query_row['stop']
                )
        db_selstr = "(resi {}-{} and chain {})".format(
                start, stop,
                df_row['chain']
                )
        db_sels.append(db_selstr)

    transform = numeric.Transformation(df_vectors, query_vectors)

    db_selstr = 'db and ('+ db_sels[0]
    for i in range(1, len(db_sels)):
        db_selstr += ' or '
        db_selstr += db_sels[i]
    db_selstr += ')'

    query_selstr = query_selstr[:-4]
    query_selstr += ' and chain A'

    print(query_selstr)
    print(db_selstr)

    subprocess.call(['./launch_pymol.sho', query_path + '/', 
        pdb, query_selstr,
        db_selstr])


def test():
    args = docopt.docopt(__doc__)
    results = pd.read_pickle(args['<results>'])
    df = pd.read_pickle(args['--dataframe'])
    # Re-build query dataframe
    init()
    helices = pd.read_pickle('test_files/boundary/cluster_representatives/query_helices.pkl')

    results = results.sort_values(by='matches')
    testrow = results.iloc[-5]
    testgraph = testrow['graph']

    session_from_graph(testgraph, helices, df)

if __name__=='__main__':
    test()
