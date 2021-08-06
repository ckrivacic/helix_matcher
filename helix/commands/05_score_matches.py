'''
Usage:
    helix 05_score_matches <workspace> [options]

options:
    --dataframe=PKL, -d
        Path to dataframe of helix vectors  [default: nr_dataframes/final.pkl]
'''
from helix.analysis import clash
from helix import workspace as ws
import docopt
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

def session_from_graph(workspace, results_row, query_df, db_df, alpha):

    def chain_from_name(string):
        chainno = int(string.split('_')[-1])
        chain = chr(ord('@')+chainno)
        return chain

    clash_score = clash.ClashScore(workspace, results_row, db_df, query_df,
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


def match_scaffold(workspace, filepath):
    '''Take a filename and find out which helix scaffold it belongs to'''
    filename = os.path.basename(filepath)
    for helix in workspace.scaffolds:
        basename = workspace.basename(helix)
        if filename.startswith(basename):
            return basename
    return None


def score_matches(workspace, results, query_df, db_df):
    '''Go through results dataframe and score the matches'''
    # for i in range(0, 100): # Review top 100 matches for now.
        # testrow = results.iloc[i]
    alphapath = query_df.iloc[0]['path']
    alphapath = clash.get_relative_path(workspace, alphapath)
    alpha = clash.get_alphashape(alphapath, plot=True)
    query_xyz = workspace.query_CAs
    results['clash_score'] = None
    results['rifdock_score'] = None
    numrows = results.shape[0]
    curr = 0

    for idx, row in results.iterrows():
        curr += 1
        print('Row {} out of {}'.format(curr, numrows))
        clash_score = clash.ClashScore(workspace, row, db_df, query_df,
                alpha=alpha)
        clash_score.apply()
        print('CLASH SCORE IS {}'.format(clash_score.score))
        results.at[idx,'clash_score'] = clash_score.score
        rifdock_score = 0
        if clash_score.subgraph:
            for node in clash_score.subgraph:
                query_idx = node[1]
                query_row = query_df.loc[query_idx]
                print(query_row)
                helixpath = os.path.realpath(query_row['path'])
                helixpath = clash.get_relative_path(workspace, helixpath)
                print('-----------------')
                # turnno = os.path.basename(helixpath).split('_')[0][0]
                score_file = \
                        match_scaffold(workspace, helixpath)\
                        + '.scores'
                scorepath = os.path.join(
                        os.path.dirname(query_row['path']),
                        score_file
                        )
                scorepath = clash.get_relative_path(workspace, scorepath)
                print(scorepath)
                # In the future, would be more efficient to just open all
                # the score files once and save them to a dataframe.
                with open(scorepath, 'r') as f:
                    for line in f:
                        line = line.strip('\n')
                        if line.startswith(os.path.basename(helixpath)):
                            print(line)
                            score_line = line
                            break
                score = float(score_line.split()[11])
                rifdock_score += score
            print('RIFDOCK SCORE IS {}'.format(rifdock_score))
            print('====================================')
            results.at[idx,'rifdock_score'] = rifdock_score
        # print(row)
        # print(results.iloc[idx])

    return results


def test_scoring():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    for target in workspace.targets:
        match_workspace = \
                ws.workspace_from_dir(workspace.target_match_path(target))
        for result in match_workspace.outputs:
            try:
                output = pd.read_pickle(result)
            except:
                with open(result, 'rb') as f:
                    output = pickle.load(f)

            suffix = result.split('.')[0].split('_')[-1]
            # output = pd.read_pickle('matcher_outputs/query_results_000.pkl')
            # helixpath = os.path.join(
                        # os.path.expanduser('~/intelligent_design/helix_matcher')
                        # )
            # helices = pd.read_pickle(
                    # os.path.join(helixpath,
                        # 'rifdock/boundary/cluster_representatives/4_turn/query_helices.pkl')
                    # )
            try:
                helices = pd.read_pickle(match_workspace.all_scaffold_dataframe)
            except:
                with open(match_workspace.all_scaffold_dataframe, 'rb') as f:
                    helices = pickle.load(f)
            # df = pd.read_pickle('dataframes_clash/final.pkl')
            df = pd.read_pickle(match_workspace.dataframe_path)
            results = score_matches(match_workspace, output, helices, df)
 
            results.to_pickle('results_scored_{}.pkl'.format(suffix))


def test():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    # results = pd.read_pickle(args['<results>'])
    # df = pd.read_pickle(args['--dataframe'])
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

        session_from_graph(workspace, testrow, helices, df, alpha)

def main():
    test_scoring()

if __name__=='__main__':
    # test()
    test_scoring()
