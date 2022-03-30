'''
Make a PyMOL session for each result, one at a time, in ascending score
order.

Usage:
    helix export_matches <workspace> [options]

Options:
    --target=STR, -t  Only view matches for the given target
'''
from pyrosetta.rosetta.core.pose import append_pose_to_pose
from helix import workspace as ws
import docopt
import pandas as pd
import pickle5 as pickle
import pymol
from pyrosetta import pose_from_file
from pyrosetta import init
from helix.utils import numeric
from helix.utils import utils
import os


def get_pymol_transform(transformation):
    pymol_transform = []

    for i in range(0, 3):
        for item in transformation.rotation[i]:
            pymol_transform.append(str(item))
        pymol_transform.append(str(transformation.translation[i]))
    # pymol_transform.extend([0,0,0,0])
    
    return pymol_transform


def session_from_graph(match_workspace, results_row, query_df, db_df):

    init()

    subgraph = results_row['subgraph']

    df_row = db_df.loc[subgraph[0][0]]
    if 'path' in df_row:
        pdb = df_row['path']
        dfpose = pose_from_file(pdb)
    else:
        pdb = df_row['name'].split('_')[0]
        pdb_name = pdb
        pymol.cmd.fetch(pdb, pdb_name)
        dfpose = pose_from_file(pdb + '.cif')
    query_vectors = []
    df_vectors = []

    match_results = ['rmsd', 'clash_score', 'interweave_score', 'total_match_score']
    # Put files in os.path.join(workspace.focus_dir, 'pdbs')
    pdbdir = os.path.join(match_workspace.focus_dir, 'pdbs')
    os.makedirs(pdbdir, exist_ok=True)
    i = 0
    filename = df_row['name']
    filepath = os.path.join(pdbdir, filename + f'_{i}.pdb.gz')
    while os.path.exists(filepath):
        i += 1
        filepath = os.path.join(pdbdir, filename + f'_{i}.pdb.gz')

    query_rows = []
    for node in subgraph:
        df_idx = node[0]
        query_idx = node[1]
        df_row = db_df.loc[df_idx]
        query_row = query_df.loc[query_idx]
        for result_type in match_results:
            query_row[result_type] = results_row[result_type]
        # To get (hopefully) accurate accounting of residues before & after combining poses, find the residue
        # index for the specific chain and use that, since we will always make the new scaffold chain A
        # Note that we only want to do this if split_chains was not selected during the scan_helices step,
        # so we have to implement a check to see if the residues are part of the correct chain or not.
        # Also note that this is intended to be used with split_chains ON, so in principle this shouldn't be
        # necessary; however, in the future, it could be useful (or at least interesting) to match multi-chain units,
        # and this avoids one potential problem with that (though many other problems will still have to be tackled).
        # The intention is for these 'rosetta_helix_start/stop' values to be used to transfer helix restypes over from
        # the designed docked helices.
        rosetta_chain = dfpose.chain(df_row['start'])
        pdb_chain = df_row['chain']
        if ' ' in pdb_chain:
            pdb_chain = pdb_chain.split(' ')[1]
        if dfpose.pdb_info().pose2pdb(dfpose.chain_begin(rosetta_chain)).split(' ')[1] == pdb_chain:
            query_row['rosetta_helix_start'] = df_row['start'] - dfpose.chain_begin(rosetta_chain) + 1
            query_row['rosetta_helix_stop'] = df_row['stop'] - dfpose.chain_begin(rosetta_chain) + 1
        else:
            query_row['rosetta_helix_start'] = df_row['start']
            query_row['rosetta_helix_stop'] = df_row['stop']
        query_row['superimposed_file'] = os.path.relpath(os.path.abspath(filepath), workspace.root_dir)
        query_rows.append(query_row)

        query_vectors.extend(query_row['vector'])
        df_vectors.extend(df_row['vector'])

    chain = df_row['chain']
    if ' ' in chain:
        chain = chain.split(' ')[1]
    dfpose = utils.pose_get_chain(dfpose, chain)
    transform = numeric.Transformation(df_vectors, query_vectors)

    targetpose = pose_from_file(match_workspace.target_path)
    targetpose.pdb_info().set_chains('B')
    dfpose.pdb_info().set_chains('A')
    rotation = numeric.np_array_to_xyzM(transform.rotation)
    translation = numeric.np_array_to_xyzV(transform.translation)
    dfpose.apply_transform_Rx_plus_v(rotation, translation)
    append_pose_to_pose(dfpose, targetpose)

    dfpose.dump_pdb(filepath)

    return query_rows


def get_realpath(row):
    full = os.path.realpath(os.path.join(
        workspace.root_dir,
        row['path']
        )
    )
    return os.path.relpath(full, workspace.root_dir)


def main():
    args = docopt.docopt(__doc__)
    global workspace
    workspace = ws.workspace_from_dir(args['<workspace>'])

    # df is the scaffold database dataframe (NOT the individual helices)
    try:
        df = pd.read_pickle(workspace.dataframe_path)
    except:
        with open(workspace.dataframe_path, 'rb') as f:
            df = pickle.load(f)

    if args['--target']:
        targets = [args['--target']]
    else:
        targets = workspace.targets
    for target in targets:
        match_workspace = \
                ws.workspace_from_dir(workspace.target_match_path(target))
        # helices is the dataframe of docked helix orientations
        try:
            helices = pd.read_pickle(match_workspace.all_scaffold_dataframe)
        except:
            with open(match_workspace.all_scaffold_dataframe, 'rb') as f:
                helices = pickle.load(f)
        scored_outputs = match_workspace.scored_outputs
        results = pd.DataFrame()
        for outfile in scored_outputs:
            try:
                scores = pd.read_pickle(outfile)
            except:
                with open(outfile, 'rb') as f:
                    scores = pickle.load(f)
            results = pd.concat([results, scores],
                    ignore_index=True)
        bins = [0, 10, 20, 30, 40, 50]
        results['binned_match_score'] = \
            pd.cut(results['total_match_score'], bins)
        results = results[results.name != '3g67_1']
        results = results.sort_values(by=['n_matched_helices',
                                          'binned_match_score',
                                          'interface_score_sum',
                                          'rmsd'], ascending=True)

        # Get designed helix scores and merge with results
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)
        designed_helix_scores = rif_workspace.get_scores()
        if 'design_file' not in designed_helix_scores.columns:
            designed_helix_scores['design_file'] = designed_helix_scores['patchman_file']

        helices.drop('name', axis=1, inplace=True)
        helices['design_file'] = helices.apply(get_realpath, axis=1)
        helices['copy_index'] = helices.index
        helices = helices.merge(designed_helix_scores, how='left', on='design_file').set_index('copy_index')

        exported_df = []
        for i in range(0, 100):
            if i < results.shape[0]:
                testrow = results.iloc[i]
                print(testrow)
                exported_df.extend(session_from_graph(match_workspace, testrow, helices, df))

        final = pd.DataFrame(exported_df)
        final.set_index('design_file')
        final.to_pickle(os.path.join(match_workspace.focus_dir, 'pdbs', 'exported_pdbs.pkl'))

