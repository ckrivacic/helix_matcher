'''
Usage:
    score_matches <workspace> [options]

options:
    --dataframe=PKL, -d
        Path to dataframe of helix vectors  [default: nr_dataframes/final.pkl]
    --plot-alphashape, -p
        Show a plot that displays the alphashape of the target protein
    --task=INT, -t
        Specify a task number if running locally  [default: 0]
    --ntasks=INT, -n
        How many tasks to split each result dataframe into  [default: 1]
    --target=TAR, -t
        Only run scoring analysis for a certain target
    --threshold=NUM
        Threshold over which results will not be saved.  [default: 50]
'''
from helix.analysis import clash
from helix.commands.view_matches import get_pymol_transform
from helix import workspace as ws
import docopt
import prody
import pandas as pd
import pickle5 as pickle
import sys
import networkx as nx
from pyrosetta import pose_from_file
from pyrosetta import init
from helix.utils import numeric
import subprocess
import os

import traceback


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

    init()
    def chain_from_name(string):
        chainno = int(string.split('_')[-1])
        chain = chr(ord('@')+chainno)
        return chain

    clash_score = clash.Score(workspace, results_row, db_df, query_df,
            alpha=alpha)
    clash_score.apply()
    print('SCORE IS {}'.format(clash_score.score), flush=True)
    if clash_score.score > 150:
        return
    subgraph = clash_score.subgraph
    print(subgraph, flush=True)
    query_selstr = ""
    db_selstr = "db and "

    query_path = os.path.dirname(query_df.loc[0]['path'])

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

    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]

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
    print(db_selstr, flush=True)

    query_selstr = query_selstr[:-4]
    query_selstr += ' and chain A'

    qobjs = qobjs[:-3] + ')'

    print(query_selstr, flush=True)
    print(db_selstr, flush=True)

    script_path = os.path.dirname(os.path.realpath(__file__))
    script_path = os.path.join(script_path, '..', 'analysis',
            'launch_pymol.sho')
    cmd = [script_path]
    cmd.append(query_path + '/')
    cmd.append(pdb)
    cmd.append(query_selstr)
    cmd.append(db_selstr)
    cmd.append(df_row['chain'])
    cmd.extend(pymol_transform)
    cmd.append(qobjs)

    print(cmd, flush=True)

    subprocess.call(cmd)


def collect_scores(results_row, query_df):
    print(results_row, flush=True)
    query_path = os.readlink(query_df.loc[0]['path'])


def match_scaffold(workspace, filepath):
    '''Take a filename and find out which helix scaffold it belongs to'''
    filename = os.path.basename(filepath)
    for helix in workspace.scaffolds:
        basename = workspace.basename(helix)
        if filename.startswith(basename):
            return basename
    return None


def apply(scorer, cutoff=50):
    '''
    Find the score for each subgraph, return the score and subgraph
    of the best-scoring transformation
    '''
    best_score = 9999
    best_interweave_score = 9999
    best_clash_score = 9999
    best_subgraph = None
    print(f"PDB PATH: {scorer.pdb_path}", flush=True)
    print(f"CHAIN: {scorer.chain}")
    original_atoms = prody.parsePDB(scorer.pdb_path,
            chain=scorer.chain).select('backbone')
    subgraphs = []
    rows = []
    for subgraph in scorer.subgraphs:
        atoms = original_atoms.copy()
        # df_vectors, query_vectors = self.get_vectors(subgraph)
        df_rows, query_rows = scorer.get_helix_rows(subgraph)
        df_vectors = scorer.get_vectors(df_rows)
        query_vectors = scorer.get_vectors(query_rows)
        transform = numeric.Transformation(df_vectors, query_vectors)
        print(f'ROTATION: {transform.rotation}', flush=True)
        print(f'TRANSLATION: {transform.translation}', flush=True)
        print(f'PYMOL TRANFORMATION: {get_pymol_transform(transform)}', flush=True)

        # Don't like that I have to use a try/except here, but I'm having some trouble in a test case, potentially
        # due to it being an old database? Shouldn't see these errors when matching to a de novo scaffold library.
        try:
            prody_transform =\
                    prody.measure.transform.Transformation(transform.rotation,
                    transform.translation)
            prody_transform.apply(atoms)
            score = scorer.calculate(atoms)
            interweave_score = scorer.calc_interweave_score(atoms, df_rows, query_rows)
            rmsd = scorer.calc_rmsd(atoms, df_rows, query_rows)
            parallel_rmsd, rmsd_lengths = scorer.parallel_rmsd(atoms, df_rows, query_rows, transform)
        except Exception as e:
            print("\033[91mError ocurred for the following.")
            print(f"\033[91mPDB: {scorer.pdb_path}")
            print(f"\033[91mChain: {scorer.chain}")
            print("\033[91mDF rows:")
            print('\033[0m')
            print(df_rows)
            print("\033[91mQuery rows:")
            print('\033[0m')
            print(query_rows)
            print('Error was:')
            print(f"\033[91m{e}\033[0m")
            print(traceback.format_exc())
            continue
        # assert(self.interweave_score is not None)

        ## print(self.interweave_score)
        print(f'Current clash score: {score}', flush=True)
        print(f'Current RMSD: {rmsd}', flush=True)
        print(f'Current parallel RMSD: {parallel_rmsd}', flush=True)

        if score + interweave_score < best_score:
            best_interweave_score = interweave_score
            best_clash_score = score
            best_score = score + interweave_score
            best_subgraph = subgraph
        if score + interweave_score < cutoff:
            subgraphs.append(subgraph)
            row = {
                    'name': scorer.name,
                    'path': scorer.pdb_path,
                    'chain': scorer.chain,
                    'clash_score': score,
                    'interweave_score': interweave_score,
                    'rmsd': rmsd,
                    'parallel_rmsd': parallel_rmsd,
                    'parallel_rmsd_lengths': rmsd_lengths,
                    'total_match_score': score + interweave_score,
                    'subgraph': subgraph,
                    'n_matched_helices': len(df_rows),
                    }
            # rosetta_score = 0
            tracked_scores = ['interface_score', 'n_hbonds', 'size', 'shape_complementarity',
                              'contact_molecular_surface', 'buns_all', 'buns_sc', 'buried_npsa_helix', 'delta_buried_npsa',
                              'exposed_hydrophobics', 'packstat', 'percent_helical', 'n_interface_residues', 'complex_sasa',
                              'delta_sasa', 'crossterm_energy', 'interface_packstat', 'delta_unsat', 'interface_dG',
                              'sequence_identity',]
            scoredict = {}
            for sc in tracked_scores:
                scoredict[sc] = 0

            # length = 0
            # rifdock_score = 0
            for node in subgraph:
                query_idx = node[1]
                query_row = scorer.query_helices.loc[query_idx]
                # Helixpath will follow symlink.
                helixpath = os.path.realpath(os.path.join(scorer.workspace.root_dir, query_row['path']))
                print(f'helixpath: {helixpath}', flush=True)
                # helixpath = clash.get_relative_path(workspace, helixpath)
                # turnno = os.path.basename(helixpath).split('_')[0][0]
                '''
                # RIFDock stuff
                score_file = \
                        match_scaffold(scorer.workspace, helixpath)\
                        + '.scores'
                scorepath = os.path.join(
                        os.path.dirname(query_row['path']),
                        score_file
                        )
                scorepath = clash.get_relative_path(scorer.workspace, scorepath)
                '''
                # rosetta_scorepath = os.path.join(
                #         os.path.dirname(query_row['path']),
                #         'scores.pkl'
                #         )
                # rosetta_scorepath = os.path.join(
                #     scorer.workspace.cluster_outputs, '..', 'scores','final.pkl'
                # )
                # try:
                #     rosetta_scores = pd.read_pickle(rosetta_scorepath)
                # except:
                #     with open(rosetta_scorepath, 'rb') as f:
                #         rosetta_scores = pickle.load(f)
                rosetta_scores = scorer.workspace.rifdock_workspace.get_scores()
                # In the future, would be more efficient to just open all
                # the score files once and save them to a dataframe.
                '''
                # More RIFDock stuff
                with open(scorepath, 'r') as f:
                    for line in f:
                        line = line.strip('\n')
                        if line.startswith(os.path.basename(query_row['path'])):
                            score_line = line
                            break
                rifdock_score += float(score_line.split()[11])
                '''
                # print(helixpath)
                # print(rosetta_scores)
                if 'design_file' not in rosetta_scores.columns:
                    rosetta_scores['design_file']  = rosetta_scores['patchman_file']
                # I should eventually change it so that 'design_file' is relative to workspace root directory,
                # as this would be safer than looking for the patchman file in case one wants ot output multiple designs per helix.
                score_row =\
                        rosetta_scores[rosetta_scores['patchman_file']==os.path.relpath(
                                helixpath,
                                scorer.workspace.root_dir)]
                for sc in scoredict:
                    scoredict[sc] += score_row.iloc[0][sc]
                # rosetta_score += float(score_row['interface_score'])
                # length += float(score_row['size'])
            for sc in scoredict:
                row[sc + '_sum'] = scoredict[sc]
            # row['rosetta_score'] = rosetta_score
            # row['rifdock_score'] = rifdock_score
            # row['per_residue_interface_score'] = rosetta_score / length
            # print('RIFDOCK SCORE IS {}'.format(rifdock_score))
            # print('ROSETTA SCORE IS {}'.format(rosetta_score))
            # # print('ROSETTA NORMALIZED SCORE IS {}'.format(rosetta_score
            #     / length))
            print('CLASH SCORE IS {}'.format(score), flush=True)
            print('INTERWEAVE SCORE IS {}'.format(interweave_score), flush=True)
            # session_from_graph(workspace, row, query_df, db_df, alpha)
            print('====================================', flush=True)
            rows.append(row)

    scorer.interweave_score = best_interweave_score
    scorer.score = best_clash_score
    scorer.subgraph = best_subgraph
    return rows


def score_matches(workspace, results, query_df, db_df, plot=False,
        threshold=50, verbose=False):
    '''Go through results dataframe and score the matches'''
    # for i in range(0, 100): # Review top 100 matches for now.
        # testrow = results.iloc[i]
    alphapath = query_df.iloc[0]['path']
    alphapath = clash.get_relative_path(workspace, alphapath)
    alpha = clash.get_alphashape(alphapath, plot=plot)
    # query_xyz = workspace.query_CAs
    results['clash_score'] = None
    results['rifdock_score'] = None
    numrows = results.shape[0]
    curr = 0

    scored_results = []
    for idx, row in results.iterrows():
        curr += 1
        print('Row {} out of {}'.format(curr, numrows), flush=True)
        clash_score = clash.Score(workspace, row, db_df, query_df,
                alpha=alpha)
        # clash_score.apply()
        candidates = apply(clash_score, cutoff=threshold)
        if verbose:
            print('Candidates:', flush=True)
            print(candidates, flush=True)
        scored_results.extend(candidates)

    return pd.DataFrame(scored_results)


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
                df = pd.read_pickle(match_workspace.dataframe_path)
                helices = pd.read_pickle(match_workspace.all_scaffold_dataframe)
            except:
                with open(match_workspace.all_scaffold_dataframe, 'rb') as f:
                    helices = pickle.load(f)
                with open(match_workspace.dataframe_path, 'rb') as f:
                    df = pickle.load(f)
            # df = pd.read_pickle('dataframes_clash/final.pkl')
            print(match_workspace.dataframe_path, flush=True)
            results = score_matches(match_workspace, output, helices,
                    df, plot=args['--plot-alphashape'])
 
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
    # test_scoring()
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    if not hasattr(workspace, 'relative_orientations'):
        print('ERROR: Not a match workspace.', flush=True)
        sys.exit()
    # targets = workspace.targets
    if 'SGE_TASK_ID' in os.environ:
        task = int(os.environ['SGE_TASK_ID']) - 1
    else:
        task = int(args['--task']) - 1
    ntasks = int(args['--ntasks'])
    # if args['--target']:
    #     match_workspace = \
    #             ws.workspace_from_dir(workspace.target_match_path(args['--target']))
    #     dataframes = match_workspace.outputs
    # else:
    dataframes = workspace.outputs
    num_dataframes = len(dataframes)
    this_job = task % num_dataframes
    subjob = task // num_dataframes

    # match_workspace = ws.workspace_from_dir(os.path.dirname(dataframes[this_job]))
    result = dataframes[this_job]
    print(f'Loading dataframe {result}', flush=True)

    try:
        output = pd.read_pickle(result)
    except:
        with open(result, 'rb') as f:
            output = pickle.load(f)

    output_size = output.shape[0]
    interval = output_size // ntasks + 1
    # print(task)
    # print(num_dataframes)
    # print(subjob)
    start = subjob * interval
    stop = start + interval
    print("START ROW: {}".format(start), flush=True)
    print("STOP ROW: {}".format(stop), flush=True)
    output = output.iloc[start:stop]

    suffix = '_'.join(result.split('.')[0].split('_')[-2:])

    try:
        df = pd.read_pickle(workspace.dataframe_path)
        helices = pd.read_pickle(workspace.all_scaffold_dataframe)
    except:
        with open(workspace.dataframe_path, 'rb') as f:
            df = pickle.load(f)
        with open(workspace.all_scaffold_dataframe, 'rb') as f:
            helices = pickle.load(f)
    print('Using the following helix dataframe: {}'.format(
        workspace.dataframe_path), flush=True)

    results = score_matches(workspace, output, helices, df,
            plot=args['--plot-alphashape'],
            threshold=float(args['--threshold']))
    out = os.path.join(workspace.output_dir,'results_scored_{}_{}.pkl'.format(suffix, subjob))
    print('SAVING RESULTS TO {}'.format(out), flush=True)
    results.to_pickle(out)


if __name__=='__main__':
    # test()
    # test_scoring()
    main()
