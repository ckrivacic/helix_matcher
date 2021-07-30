'''
Usage:
    review_matches.py <results> [options]

options:
    --dataframe=PKL, -d
        Path to dataframe of helix vectors  [default: nr_dataframes/final.pkl]
'''
import clash
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

def get_pymol_transform(transformation):
    pymol_transform = []

    for i in range(0, 3):
        for item in transformation.rotation[i]:
            pymol_transform.append(str(item))
        pymol_transform.append(str(transformation.translation[i]))
    # pymol_transform.extend([0,0,0,0])
    
    return pymol_transform

def session_from_graph(results_row, query_df, db_df, alpha):

    def chain_from_name(string):
        chainno = int(string.split('_')[-1])
        chain = chr(ord('@')+chainno)
        return chain

    clash_score = clash.ClashScore(results_row, db_df, query_df,
            alpha=alpha)
    clash_score.apply()
    print('SCORE IS {}'.format(clash_score.score))
    if clash_score.score > 150:
        return
    subgraph = clash_score.subgraph
    print(subgraph)
    query_selstr = ""
    db_selstr = "db and "

    query_path = os.path.dirname(query_df.loc[0]['path'])

    db_sels = []
    df_row = db_df.loc[subgraph[0][0]]
    pdb = df_row['name'].split('_')[0]
    pymol.cmd.fetch(pdb)
    dfpose = pose_from_file(pdb + '.cif')
    query_vectors = []
    df_vectors = []
    qobjs = '('

    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]

        start = dfpose.pdb_info().pose2pdb(df_row['start']).split(' ')[0]
        stop = dfpose.pdb_info().pose2pdb(df_row['stop']).split(' ')[0]
        # start = df_row['start']
        # stop = df_row['stop']
        
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

    cmd = ['./launch_pymol.sho']
    cmd.append(query_path + '/')
    cmd.append(pdb)
    cmd.append(query_selstr)
    cmd.append(db_selstr)
    cmd.append(df_row['chain'])
    cmd.extend(pymol_transform)
    cmd.append(qobjs)

    print(cmd)

    subprocess.call(cmd)


def collect_scores(results_row, query_df):
    print(results_row)
    query_path = os.readlink(query_df.loc[0]['path'])


def score_matches(results, query_df, db_df):
    '''Go through results dataframe and score the matches'''
    # for i in range(0, 100): # Review top 100 matches for now.
        # testrow = results.iloc[i]
    alphapath = query_df.iloc[0]['path']
    alpha = clash.get_alphashape(alphapath)
    results['clash_score'] = None
    results['rifdock_score'] = None
    numrows = results.shape[0]
    curr = 0

    for idx, row in results.iterrows():
        curr += 1
        print('Row {} out of {}'.format(curr, numrows))
        clash_score = clash.ClashScore(row, db_df, query_df,
                alpha=alpha)
        clash_score.apply()
        print('CLASH SCORE IS {}'.format(clash_score.score))
        results.at[idx,'clash_score'] = clash_score.score
        rifdock_score = 0
        if clash_score.subgraph:
            for node in clash_score.subgraph:
                query_idx = node[1]
                query_row = query_df.loc[query_idx]
                helixpath = os.path.realpath(query_row['path'])
                turnno = os.path.basename(helixpath).split('_')[0][0]
                scorepath = os.path.join(
                        os.path.dirname(query_row['path']),
                            '{}turn.scores'.format(turnno)
                        )
                # In the future, would be more efficient to just open all
                # the score files once and save them to a dataframe.
                with open(scorepath, 'r') as f:
                    for line in f:
                        line = line.strip('\n')
                        if line.endswith(os.path.basename(helixpath)):
                            score_line = line
                            break
                score = float(score_line.split()[10])
                rifdock_score += score
            print('RIFDOCK SCORE IS {}'.format(rifdock_score))
            results.at[idx,'rifdock_score'] = rifdock_score
        # print(row)
        # print(results.iloc[idx])

    return results


def test_scoring():
    results = pd.read_pickle('matcher_outputs/query_results_000.pkl')
    helixpath = os.path.join(
                os.path.expanduser('~/intelligent_design/helix_matcher')
                )
    helices = pd.read_pickle(
            os.path.join(helixpath,
                'rifdock/boundary/cluster_representatives/4_turn/query_helices.pkl')
            )
    df = pd.read_pickle('dataframes_clash/final.pkl')
    results = score_matches(results, helices, df)
    
    results.to_pickle('results_scored.pkl')


def test():
    args = docopt.docopt(__doc__)
    results = pd.read_pickle(args['<results>'])
    df = pd.read_pickle(args['--dataframe'])
    # Re-build query dataframe
    init()
    helixpath = os.path.expanduser('~/intelligent_design/helix_matcher')
    helices = pd.read_pickle(os.path.join(helixpath,
        'test_files/boundary/cluster_representatives/query_helices.pkl'))

    results = results.sort_values(by='matches', ascending=False)
    alphapath = helices.iloc[0]['path']
    alpha = clash.get_alphashape(alphapath, chain='B')
    for i in range(0, 100): # Review top 100 matches for now.
        testrow = results.iloc[i]
        # testgraph = testrow['graph']

        session_from_graph(testrow, helices, df, alpha)

if __name__=='__main__':
    # test()
    test_scoring()
