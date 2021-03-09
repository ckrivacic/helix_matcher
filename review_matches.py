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
    query_selstr = "query and ("
    db_selstr = "db and ("

    print(query_df)

    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]

        query_selstr += "(resi {}-{} and chain {}) or ".format(
                query_row['start'], query_row['stop'], query_row['chain']
                )
        db_selstr += "(resi {}-{} and chain {}) or ".format(
                df_row['start'], df_row['stop'],
                chain_from_name(df_row['name'])
        )

    query_selstr += ")"
    db_selstr += ")"

    print(query_selstr)
    print(db_selstr)

    pymol.finish_launching()



def test():
    args = docopt.docopt(__doc__)
    results = pd.read_pickle(args['<results>'])
    df = pd.read_pickle(args['--dataframe'])
    # Re-build query dataframe
    init()
    pose = pose_from_file('test_files/boundary/cluster_representatives/combined.pdb')
    scanner = scan_helices.PoseScanner(pose)
    helices = scanner.scan_pose_helices(name='query',
            split_chains=False)
    helices = pd.DataFrame(helices)

    testrow = results.iloc[0]
    testgraph = testrow['graph']

    session_from_graph(testgraph, helices, df)


if __name__=='__main__':
    test()
