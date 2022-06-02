'''
Filter by interface score.

Usage:
    helix 04_filter_designs <workspace> <output_dir> [options]

Options:
    --filter=FILENAME, -f  A yaml file describing filters. If not provided, a
    default set of filters will be used.
    --target=STR, -t  Only run filters for a specific target
    --clear, -o  Delete existing filtered symlinks prior to running
    --copy  Copy files instead of symlinking them (NOT IMPLEMENTED)
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
        match_workspace = ws.MatchWorkspace(workspace.root_dir, target)

        if args['--clear']:
            match_workspace.clear_cluster_outputs()

        # Import all "score" files
        scores = match_workspace.get_scores()

        # Filter the final dataframe
        if args['--filter']:
            with open(args['--filter'], 'r') as f:
                filters = yaml.load(f.read())
        else:
            filters = {
                    'threshold': {
                        'percent_helical': ['>', 0.7],
                        'buns_all': ['<', 2],
                        },
                    'percentile': {
                        'interface_dG': ['<', 0.4]
                        'contact_molecular_surface': ['>', 0.5]
                        }
                    }
        print('FILTERS:')
        print(filters)
        scorelist = []
        # for name, group in scores.groupby(['patch_len']):
        #     scorelist.append(parse_filter(filters, group))
        scores = parse_filter(filters, scores)

        # Symlink 
        # outdir = match_workspace.cluster_outputs
        outdir = args['<output_dir>']
        if not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        for idx, row in scores.iterrows():
            # patchman_length = os.path.dirname(row['design_file']).split('/')[-2]
            # symoutdir = os.path.join(outdir, patchman_len)
            symlink_path = os.path.join(outdir, os.path.basename(row['design_file']))
            if os.path.islink(symlink_path):
                os.remove(symlink_path)
            relpath = os.path.relpath(os.path.join(match_workspace.root_dir, row['design_file']), outdir)
            os.symlink(relpath, symlink_path)