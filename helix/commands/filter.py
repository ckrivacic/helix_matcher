'''
Filter by interface score.

Usage:
    helix filter <workspace> [options]

Options:
    --filter, -f  A yaml file describing filters. If not provided, a
    default set of filters will be used.
    --target=STR, -t  Only run filters for a specific target
    --clear, -o  Delete existing filtered symlinks prior to running
'''
import docopt
import yaml
from helix import workspace as ws
import pandas as pd
import os
import copy


def parse_filter(filter_dict, df):
    '''Parse filters. Filters should have the following dictionary format:
        
        {
        'threshold':
            {'n_hbond': ['>', 2],}, # n_hbond should be greater than 2
        'percentile':
            {'interface_score_y': ['<', 0.5],} # Only keep data below
            # the 50th percentile in interface_score_y
        }
    '''
    orig_df = copy.deepcopy(df)
    if 'threshold' in filter_dict:
        for thresh in filter_dict['threshold']:
            operator = filter_dict['threshold'][thresh][0]
            if operator == '>':
                df = df[df[thresh] > filter_dict['threshold'][thresh][1]]
            elif operator == '<':
                df = df[df[thresh] < filter_dict['threshold'][thresh][1]]
            else:
                print('Filter dictionary had unexpected format. Please only use \'>\' or \'<\' for thresholds')

    if 'percentile' in filter_dict:
        for percentile in filter_dict['percentile']:
            operator = filter_dict['percentile'][percentile][0]
            if operator == '<':
                df = df[df[percentile] < orig_df[percentile].quantile(
                    filter_dict['percentile'][percentile][1]
                    )]
            else:
                df = df[df[percentile] > orig_df[percentile].quantile(
                    filter_dict['percentile'][percentile][1]
                    )]

    return df


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])

    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    # Iterate through all targets?
    for target in targets:
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)

        if args['--clear']:
            rif_workspace.clear_cluster_outputs()

        # Import all "score" files
        scores = rif_workspace.get_scores()
        if 'design_file' not in scores.columns:
            scores['design_file'] = scores['patchman_file']

        # Filter the final dataframe
        if args['--filter']:
            filters = yaml.load(args['--filter'])
        else:
            filters = {
                    'threshold': {
                        'percent_helical': ['>', 0.7],
                        'n_hbonds': ['>', 1],
                        },
                    'percentile': {
                        'interface_dG': ['<', 0.5]
                        }
                    }
        scores = parse_filter(filters, scores)

        # Symlink 
        outdir = rif_workspace.cluster_outputs
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        for idx, row in scores.iterrows():
            # patchman_length = os.path.dirname(row['design_file']).split('/')[-2]
            # symoutdir = os.path.join(outdir, patchman_len)
            symlink_path = os.path.join(outdir, os.path.basename(row['design_file']))
            if os.path.islink(symlink_path):
                os.remove(symlink_path)
            relpath = os.path.relpath(os.path.join(rif_workspace.root_dir, row['design_file']), outdir)
            os.symlink(relpath, symlink_path)