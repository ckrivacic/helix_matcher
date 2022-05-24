'''
Filter by interface score.

Usage:
    helix filter <workspace> [options]

Options:
    --filter=FILENAME, -f  A yaml file describing filters. If not provided, a
    default set of filters will be used.
    --target=STR, -t  Only run filters for a specific target
    --clear, -o  Delete existing filtered symlinks prior to running
'''
import docopt
import yaml
from helix import workspace as ws
import pandas as pd
import os

from helix.utils.utils import parse_filter


def get_patch_length(row):
    return row['patchman_file'].split('/')[2]


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])

    if args['--target']:
        targets = [workspace.target_rifdock_path(args['--target'])]
    else:
        targets = workspace.all_rifdock_workspaces

    # Iterate through all targets?
    for target in targets:
        print(f'Filtering target {target}')
        rif_workspace = ws.RIFWorkspace(workspace.root_dir, target)

        if args['--clear']:
            rif_workspace.clear_cluster_outputs()

        # Import all "score" files
        scores = rif_workspace.get_scores()
        if 'patch_len' not in scores.columns:
            scores['patch_len'] = scores.apply(get_patch_length, axis=1)
        if 'design_file' not in scores.columns:
            scores['design_file'] = scores['patchman_file']

        # Filter the final dataframe
        if args['--filter']:
            with open(args['--filter'], 'r') as f:
                filters = yaml.load(f.read())
        else:
            filters = {
                    'threshold': {
                        'percent_helical': ['>', 0.7],
                        'n_hbonds': ['>', 1],
                        },
                    'percentile': {
                        'interface_dG': ['<', -0.5]
                        }
                    }
        print('FILTERS:')
        print(filters)
        scorelist = []
        for name, group in scores.groupby(['patch_len']):
            scorelist.append(parse_filter(filters, group))
        scores = pd.concat(scorelist)

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