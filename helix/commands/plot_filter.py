'''
Usage:
    helix plot_filter <workspace> [options]

Options:
    --filters=PATH, -f  Yaml file with filters  [default: project_params/filters.json]
    --trim=LENGTH  Trim benchmark to only have helices of above a certain length
'''

import docopt
import os
import yaml
from helix.utils import utils
import helix.workspace as ws
from helix.commands.filter import parse_filter


def load_filters(path):
    with open(path, 'r') as f:
        yaml_filters = yaml.load(f.read())
    return yaml_filters


def get_coverage(patch, benchmark, args):
    recapitulated_benchmark_data = get_benchmark_results(patch, benchmark, args)
    return recapitulated_benchmark_data[recapitulated_benchmark_data['rmsd'] < 1].shape[0] / \
                recapitulated_benchmark_data.shape[0]


def get_accuracy(patch, filter_name):
    outrows = []
    for group_name, group in patch.groupby(['name', 'target']):
        name, target = group_name
        fraction_subA = group[group['best_rmsd'] < 1].shape[0] / group.shape[0]
        outrows.append({
            'name': name,
            'target': target,
            'fraction_subA': fraction_subA
        })

    dataframe = pd.DataFrame(outrows)
    dataframe['filter'] = filter_name



def get_benchmark_results(patch, benchmark, args):

    '''patch = patchman design dataframe, benchmark = benchmark dataframe'''
    no_patch_match = []
    outrows = []
    if args['--length'] == '14':
        patch = patch[patch['patch_len'] == 'len_14']
    elif args['--length'] == '28':
        patch = patch[patch['patch_len'] == 'len_28']
    for name, group in benchmark.groupby(['name', 'target']):
        print(name)
        patch_group = patch[(patch['name'] == name[0]) & (patch['target'] == name[1])]
        for idx, row in group.iterrows():
            benchmark_resis = row['start_stop']
            patch_subgroup = patch_group[patch_group['start_stop'] == benchmark_resis]
            best_patchman = patch_subgroup.sort_values(by='best_rmsd', ascending=True)
            if best_patchman.shape[0] == 0:
                no_patch_match.append((row['name'], row['target'], row['start_stop']))
            else:
                best_patchman = best_patchman.iloc[0]
                outrows.append({'name': row['name'], 'target': row.target, 'start_stop': row.start_stop,
                                'rmsd': best_patchman.best_rmsd, 'protocol': 'PatchMAN'})

    return pd.DataFrame(outrows)


def parse_filter_name(filter):
    name = ''
    types = ['threshold', 'percentile']
    for t in types:
        for f in filter[t]:
            name += f
            name += '_'
            name += filter[t][f][0]
            name += '_'
            name += str(filter[t][f][1])


def main():
    args = docopt.docopt(__doc__)
    bench_path = os.path.join(
        os.path.dirname(os.path.realpath(__file__)),
        '..', 'benchmark', 'interface_finder', 'final_consolidated.pkl'
    )

    workspace = ws.workspace_from_dir(args['<workspace>'])
    outdir = 'helper_dataframes/'
    outname = os.path.basename(args['--filters']).split('.')[0]
    if args['--trim']:
        outname += '_trimmed_{}'.format(args['--trim'])
    outname_acc = outname + '_accuracy'
    outname_cov = outname + '_coverage'
    outpath_acc = os.path.join(outdir, outname_acc + '.pkl')
    outpath_cov = os.path.join(outdir, outname_cov + '.pkl')

    if not os.path.exists(outpath_cov) or not os.path.exists(outpath_acc):
        # Load dataframe
        df_path = os.path.join(workspace.rifdock_outdir, 'combined_benchmark_rev', 'final.pkl')
        df = utils.safe_load(df_path)

        benchmark = utils.safe_load(bench_path)

        if args['--trim']:
            benchmark['helix_length'] = benchmark.apply(lambda x: len(x['pdb_resis']), axis=1)

            def count(benchdf, length=14):
                return benchdf[benchdf['helix_length'] > length].shape[0]

            trim = []
            trim_cutoff = int(args['--trim'])
            for name, group in benchmark.groupby(['name', 'target']):
                if count(group, length=trim_cutoff) < 2:
                    trim.append(name)

            print('The following targets will be trimmed:')
            print(trim)
            df['tup'] = df.apply(lambda x: (x['name'], x['target']), axis=1)
            for name in trim:
                print(df[df['tup'] == name])

            print(df.shape)
            # df = df[(~df['name'].isin(trim_name)) & (~df['target'].isin(trim_chain))]
            df = df[~df['tup'].isin(trim)]
            print(df.shape)
        filters = load_filters(args['--filters'])
        # Filter on a single metric, several cutoffs;
        filtered_dataframes = []
        filter_names = []
        for i in range(0, len(filters)):
            current_filter = filters[i]
            filtered_df = parse_filter(current_filter, df)
            filtered_dataframes.append(filtered_df)
            filter_names.append(parse_filter_name(current_filter))

        coverage_results = []
        accuracy_results = []
        coverage_results.append({
            'filter': 'None',
            'coverage': get_coverage(df),
        })
        accuracy_results.append(get_accuracy(df, 'None'))
        for filtered_df, filter_name in zip(filtered_dataframes, filter_names):
            coverage_results.append({
                'filter': filter_name,
                'coverage': get_coverage(filtered_df)
            })
            accuracy_results.append(get_accuracy(filtered_df, filter_name))
        print('COVERAGE RESULTS')
        print(coverage_results)
        print('ACCURACY RESULTS')
        print(accuracy_results)

    # For each cutoff, calculate benchmark coverage and accuracy