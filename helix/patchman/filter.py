'''
Filter by interface score.

Usage:
    helix filter <workspace> [options]

Options:
    --filter, -f  A yaml file describing filters. If not provided, a
    default set of filters will be used.
    --target, -t  Only run filters for a specific target
    --clear, -o  Delete existing filtered symlinks prior to running
'''
import docopt
import yaml
from helix import workpace as ws


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
    if 'threshold' in filter_dict:
        for thresh in filter_dict['threshold']:
            operator = filter_dict['threshold'][thresh][0]
            if operator == '>':
                df = df[df[thresh] > filter_dict['threshold'][thresh][1]]
            elif operator == '<':
                print(df['buns_all'])
                df = df[df[thresh] < filter_dict['threshold'][thresh][1]]
            else:
                print('You fucked up your filter dictionary')
                print('Please only use \'>\' or \'<\' for thresholds')

    if 'percentile' in filter_dict:
        for percentile in filter_dict['percentile']:
            operator = filter_dict['percentile'][percentile][0]
            outrows = []
            if operator == '<':
                for name, group in df.groupby(['name_x', 'chain',
                    'pdb_tup']):
                    group = group[group[percentile] < group[percentile].quantile(
                        filter_dict['percentile'][percentile][1]
                        )]
                    outrows.append(group)
            else:
                for name, group in df.groupby(['name_x', 'chain',
                    'pdb_tup']):
                    group = group[group[percentile] < group[percentile].quantile(
                        filter_dict['percentile'][percentile][1]
                        )]
                    outrows.append(group)

        return pd.concat(outrows, ignore_index=True)
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
                    'interface_score_y': ['<', 0.5]
                    }
                }
        scores = parse_filter(filters, scores)

        # Symlink 
        outpath = os.path.join(workspace.focus_dir, 'filtered')
