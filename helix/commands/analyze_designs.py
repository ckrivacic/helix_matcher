'''
Analyze different design methods.

Usage:
    helix analyze_designs <workspace> <plot_type> [options]

Options:
    --yaxis=COLUMN, -y  Which column to plot on the y-axis  
        [default: helix_percent_identity]

    --xaxis=COLUMN, -x  Which column to plot on the x-axis  
        [default: interface_score_y]

    --rmsd-cutoff=FLOAT, -c  Eliminate helices greater than this RMSD cutoff  
        [default: 0]

    --pose-score-cutoff=FLOAT, -p  Elimiate poses with a total score
        greater than this float  [default: 0]

    --focus-dir=STR, -f  Only plot for a specific focus dir or list of focus
        dirs (ex. 1b33_K,1b33_K_specialrot)

    --size=COLUMN, -s  Size scatter plot points by this column

    --hue=COLUMN, -h  Color by a column

    --id-cutoff=FLOAT  Filter out results that have a sequence identity
        within this #

    --patch-id-cutoff=FLOAT  filter out patches that have this amount of
        sequence identity

    --filter=JSON  Path to json file to filter helices on

    --buried-identity  Calculate buried sequence identity. WARNING: Do
        not use without a cutoff... slow.

    --target=STR  Only show results for a specific target

    --protocol=STR  Only show results for this protocol

    --length=#, -l  Only show results for docked helices of this length (14 orr 28)

    --perres  Plot the per-residue #

    --buns-vs-bonds  (for bar plot) plot BUNS & # HBONDS side by side

    --test  Load a small dataframe for testing

    --filter-plot  Scatterplot showing unfiltered & filtered designs, as well as metrics for benchmark interfaces
'''
import docopt
from helix.commands import filter
import os
import pandas as pd
from klab import scripting
from helix.utils import utils
from helix.utils import homology
from helix.utils import plotting
import helix.workspace as ws
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pyrosetta import *
from pyrosetta.rosetta.core.select import residue_selector
from helix.utils.colors import palette


def first3_res_correct(pose, res, seq):
    print(seq)
    print(pose.residue(res).name1())
    first = pose.residue(res).name1() == seq[0]
    second = pose.residue(res+1).name1() == seq[1]
    third = pose.residue(res+2).name1() == seq[2]
    return first and second and third


def calc_designed_identity(row):
    '''Calculate percent identity of designed residues.
    Assumes residues in row['helix_resis'] and row['design_resis'] are
    in Rosetta numbering.'''
    helix_start = row['helix_resis'][0]
    helix_stop = row['helix_resis'][0]
    seq_indices = []
    for idx, resi in enumerate(row['designed_residues']):
        if resi >= helix_start or resi <= helix_stop:
            seq_indices.append(resi - helix_start)

    seqlen = 0
    identity = 0
    for idx, aa in enumerate(row.helix_seq):
        if idx in seq_indices:
            seqlen += 1
            if row.benchmark_seq[idx] == aa:
                identity += 1

    if seqlen == 0:
        print('No designed sequences?')
        print('FILE: {}'.format(row.design_file))
        print('DESIGNED RESIS: {}'.format(row.designed_residues))
        print('HELIX RESIS: {}'.format(row.helix_resis))
        return np.nan

    return 100 * identity / seqlen


def create_web_logo(name, group, preview=False):
    '''Make a sequence logo'''
    import weblogo
    import tempfile
    sequences = group['helix_seq']
    print('NAME')
    print(name)
    print('# sequences:')
    print(len(list(sequences)))
    print('Benchmark seq:')
    print(group.iloc[0]['benchmark_seq'])
    for idx, row in group.iterrows():
        print('Residues witheld for {}'.format(row['patchman_file_x']))
        print(row['residues_witheld'])
    sequences = weblogo.seq.SeqList(
            [weblogo.seq.Seq(x) for x in sequences],
            alphabet=weblogo.seq.unambiguous_protein_alphabet,
    )

    logo_data = weblogo.LogoData.from_seqs(sequences)
    logo_options = weblogo.LogoOptions()
    logo_options.title = name
    logo_format = weblogo.LogoFormat(logo_data, logo_options)
    if preview:
        logo_file = tempfile.NamedTemporaryFile(
                'wb', prefix='weblogo_', suffix='.pdf'
                )
    else:
        outdir = os.path.join(workspace.rifdock_outdir, 'logos')
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        outname = ''
        for n in name[0:3]:
            outname += str(n)
        outname += '_{}'.format(name[3][0])
        outname += '_{}'.format(name[3][1])
        outfile = os.path.join(outdir,  outname + '.pdf')
        print('Saving to {}'.format(outfile))
        logo_file = open(outfile, 'wb')
        with open(outfile[:-len('pdf')] + 'txt', 'w') as f:
            for line in str(logo_data):
                f.write(line)

    ext = os.path.splitext(logo_file.name)[-1]
    formatters = {
            '.pdf': weblogo.pdf_formatter,
            '.svg': weblogo.svg_formatter,
            '.eps': weblogo.eps_formatter,
            '.png': weblogo.png_formatter,
            '.jpeg': weblogo.jpeg_formatter,
            '.txt': weblogo.txt_formatter,
    }
    if ext not in formatters:
        scripting.print_error_and_die("'{0}' is not a supported file format".format(ext))

    document = formatters[ext](logo_data, logo_format)
    logo_file.write(document)
    logo_file.flush()

    if preview:
        pdf = os.environ.get('PDF', 'evince'), logo_file.name
        subprocess.call(pdf)

    logo_file.close()


def calc_identity(row):
    '''Calculate sequence identity of a helix compared to the nearest
    benchmark helix'''
    seq_helix = row['helix_seq']
    seq_benchmark = row['benchmark_seq']
    assert len(seq_helix) == len(seq_benchmark)
    identity = 0
    for i, res in enumerate(seq_helix):
        if res == seq_benchmark[i]:
            identity += 1

    return 100 * (identity / len(seq_helix))


def calc_hydrophobic_identity(row):
    '''Calculate sequence identity of a helix compared to the nearest
    benchmark helix'''
    seq_helix = row['helix_seq']
    seq_benchmark = row['benchmark_seq']
    assert len(seq_helix) == len(seq_benchmark)
    identity = 0
    hydrophobics = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    n_hydro = 0
    for i, res in enumerate(seq_helix):
        if seq_benchmark[i] in hydrophobics:
            n_hydro += 1
            if res == seq_benchmark[i]:
                identity += 1

    return 100 * (identity / n_hydro)


def plot_sequence_recovery(df, args):
    '''Make a plot for sequence recovery'''
    # for col in df.columns:
        # printcol = col
        # print(df.sort_values(by=printcol)[printcol])
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    # sns.stripplot(data=df, x='protocol', y=args['--yaxis'],
            # order=order, color='.5', alpha=0.5)
    # df = df[df['protocol'] != 'deleteme']
    sns.categorical._Old_Violin = sns.categorical._ViolinPlotter
    class _My_ViolinPlotter(sns.categorical._Old_Violin):

        def __init__(self, *args, **kwargs):
            super(_My_ViolinPlotter, self).__init__(*args, **kwargs)
            self.gray = 'black'

    sns.categorical._ViolinPlotter = _My_ViolinPlotter
    means = []
    for protocol in order:
        mean = df[df['protocol']==protocol][args['--yaxis']].mean()
        means.append(mean)
    fig, ax = plt.subplots(figsize=(4,4), dpi=300)
    if args['--yaxis'] == 'n_hbonds':
        bw = 0.2
    elif args['--yaxis'] == 'helix_percent_identity':
        bw = 0.2
    elif args['--yaxis'] == 'buns_all' or args['--yaxis'] == 'buns_sc':
        bw = 0.2
    else:
        bw = None
    sns_ax = sns.violinplot(data=df, x='protocol', y=args['--yaxis'],
            order=order, bw=bw, inner=None, linewidth=1)
    plt.scatter(x=range(len(means)), y=means, c='k', marker='_', s=200)
    if args['--xaxis'] == 'protocol':
        ax.set_xticklabels(labels, ha='right')
        plt.xticks(rotation=70)
    plt.ylabel('Sequence identity (%)')
    plt.tight_layout()
    plt.show()


def multi_scatter(groups, args, ax):
    '''Plot each group separately. Up to 3. The first will have low opacity, the second normal, the third will be
    marked with Xs.'''
    colors = ["#abd0e7", palette["blue"], palette["red"]]
    i = 0
    if args['--perres']:
        x = args['--xaxis'] + '_perres'
        y = args['--yaxis'] + '_perres'
    else:
        x = args['--xaxis']
        y = args['--yaxis']
    sizes = [3, 5, 5]
    markers = ['.', 'o', 'X']
    alphas = [0.3, 1.0, 0.5]
    linewidths = [0, 1, 0]
    for group in groups:
        alpha = alphas[i]
        sns_ax = sns.scatterplot(data=group, x=x, y=y, #picker=True,
                        alpha=alpha, ax=ax, color=colors[i], s=sizes[i], marker=markers[i], #linewidth=linewidths[i]
                                 )
        i += 1
    plt.ylabel(y_names[y])
    plt.xlabel(y_names[x])

    return sns_ax


def scatterplot(df, args):
    def calc_perres(row, col):
        if type(row['patch_len']) == type('hi'):
            return row[col] / int(row['patch_len'].split('_')[1])
        else:
            return row[col] / int(row['length'])

    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    if not args['--xaxis'] == 'protocol':
        order = None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    fig, ax = plt.subplots(figsize=(4,2), dpi=300)

    if args['--perres']:
        cols = [args['--yaxis'], args['--xaxis']]
        if args['--buns-vs-bonds']:
            cols = [ 'n_hbonds', 'buns_all']
        for col in cols:
            df[f"{col}_perres"] = df.apply(calc_perres, args=(col,), axis=1)

    if args['--filter'] and args['--filter-plot']:
        import yaml
        filter_path = args['--filter']
        with open(filter_path, 'r') as f:
            filters = yaml.load(f.read())
        scorelist = []
        df_designed = df[~df['benchmark']]
        df_bench = df[df['benchmark']]
        for name, group in df_designed.groupby(['patch_len', 'name_x', 'target']):
            scorelist.append(filter.parse_filter(filters, group))
        df_filtered = pd.concat(scorelist)
        ax = multi_scatter([df_designed, df_filtered, df_bench], args, ax)

    else:
        print(df[args['--yaxis']].sort_values())
        if args['--size']:
            size = args['--size']
            ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                    hue=hue, size=size, sizes=(50,300), picker=True)
        else:
            ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                    hue=hue, picker=True)
    # click = plotting.ClickablePlot(ax, df, args, workspace)
    plt.tight_layout()

    plt.show()


def barplot_two_cols(df, y1, y2, args, ax):
    '''Plot two columns grouped'''
    xcenters = [i for i in range(len(order))]
    delt = 0.2
    xneg = [x - delt for x in xcenters]
    xpos = [x + delt for x in xcenters]
    ax1 = sns.barplot(data=df, x=args['--xaxis'], y=y1, palette=colors, order=order, ax=ax, ci=None)
    # for i, bar in enumerate(ax1.patches):
    #     bar.set_x(bar.get_x() - delt)
    ax2 = sns.barplot(data=df, x=args['--xaxis'], y=y2, palette=[palette['lightgray']], order=order, ax=ax, ci=None)
    # for i, bar in enumerate(ax2.patches):
    #     bar.set_x(bar.get_x() + delt)

    return ax1


def barplot(df, args):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    # labels = ['Base', 'BUNS penalty', 'BUNS penalty + prune', 
            # 'Freeze good rotamers', 'Special rotamer bonus']
    # if not args['--xaxis'] == 'protocol':
    #     order=None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    if args['--perres']:
        cols = [args['--yaxis']]
        if args['--buns-vs-bonds']:
            cols = [ 'n_hbonds', 'buns_all']
        for col in cols:
            df[f"{col}_perres"] = df.apply(lambda x: \
                                                       x[col] / int(x['patch_len'].split('_')[1]),
                                                       axis=1)
    fig, ax = plt.subplots(figsize=(4,3), dpi=300)
    # sns.stripplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
            # order=order, color='.5', alpha=0.5, ax=ax)
    y_axis = args['--yaxis']
    if args['--perres']:
        y_axis += '_perres'
    if args['--buns-vs-bonds']:
        # cols = ['design_file', 'protocol', 'value', 'type']
        # if args['--perres']:
        #     df_buns = df[['design_file', 'protocol', 'buns_all_perres']]
        #     df_buns['type'] = 'buns_perres'
        #     df_bonds = df[['design_file', 'protocol', 'n_hbonds_perres']]
        #     df_bonds['type'] = 'hbonds_perres'
        # else:
        #     df_buns = df[['design_file', 'protocol', 'buns_all']]
        #     df_buns['type'] = 'buns_all'
        #     df_bonds = df[['design_file', 'protocol', 'n_hbonds']]
        #     df_bonds['type'] = 'hbonds'
        # df_buns.columns = cols
        # df_bonds.columns = cols
        # df =pd.concat([df_buns, df_bonds], ignore_index=True)
        # y_axis = 'value'
        # hue = 'type'
        print(df)
        print(df.columns)
        sns_ax = barplot_two_cols(df, 'n_hbonds_perres', 'buns_all_perres', args, ax)
    else:
        sns_ax = sns.barplot(data=df, x=args['--xaxis'], y=y_axis,
                hue=hue, order=order, ax=ax, saturation=1, palette=colors)
    if args['--xaxis'] == 'protocol':
        ax.set_xticklabels(labels, ha='right')
        plt.xticks(rotation=70)
    plt.ylabel(y_names[y_axis])
    import matplotlib.patches

    for patch in sns_ax.patches:
        patch.set_edgecolor('black')
    fig.tight_layout()
    plt.xlabel(None)
    plt.show()


def protocol_vs_protocol(df, args):
    '''Plot a metric of one protocol vs another'''
    # dfx = df[df['protocol'] == args['--xaxis']]
    # dfy = df[df['protocol'] == args['--yaxis']]
    # dfx['patchman_basename'] = dfx.apply(lambda x:
            # os.path.basename(x['patchman_file']), axis=1)
    def relative_to_focusdir(row):
        patchfile = row['patchman_file_x']
        return os.path.join(*patchfile.split('/')[2:])
    print('--------------------------------------------------')
    df['patchman_basename'] = df.apply(relative_to_focusdir, axis=1)
    xname = args['--hue'] + '_{}'.format(args['--xaxis'])
    yname = args['--hue'] + '_{}'.format(args['--yaxis'])
    dfx = df[df.protocol==args['--xaxis']]
    # dfx.rename(columns={args['--hue']: xname})
    dfy = df[df.protocol==args['--yaxis']]
    # dfy.rename(columns={args['--hue']: yname})
    dfplot = pd.merge(dfx, dfy, on=['patchman_basename', 'target',
        'name_y'],
            suffixes=('_{}'.format(args['--xaxis']),
                '_{}'.format(args['--yaxis'])),
            validate='1:1')
    dfplot.to_pickle('test.pkl')

    ax = sns.scatterplot(data=dfplot, x=xname, y=yname, picker=True)
    clickable = plotting.ClickablePlot(ax, dfplot, args, workspace, pvp=True)
    minimum = min(min(dfplot[xname]), min(dfplot[yname]))
    maximum = max(max(dfplot[xname]), max(dfplot[yname]))
    plt.plot([minimum, maximum], [minimum, maximum], linewidth=2,
            color='r')
    plt.show()


def get_patchman_pdbid(row):
    return os.path.basename(row['patchman_file_x']).split('_')[1].upper()

def relative_path(path, workspace):
    pathlist = path.split('/')
    if os.path.basename(workspace.root_dir) not in pathlist:
        # If workspace isn't in pathlist, the path is already relative
        fpath = os.path.join(
                workspace.root_dir,
                patchman_file
                )
    else:
        start = pathlist.index(
                os.path.basename(workspace.root_dir)
                ) + 1
        fpath = os.path.join(
                workspace.root_dir,
                *pathlist[start:]
                )
    return fpath


def calc_buried_identity(row):

    init()
    pose = pose_from_file(relative_path(row['design_file'], workspace))
    print(relative_path(row['design_file'], workspace))
    buried = residue_selector.LayerSelector()
    buried.set_layers(True, False, False)
    buried_res = buried.apply(pose)
    buried_res = utils.res_selector_to_size_list(buried_res, pylist=True)
    start = pose.chain_begin(2)

    seq_helix = row['helix_seq']
    while not first3_res_correct(pose, start, seq_helix):
        start += 1
    # idx + start should get to the right resnum 

    seq_benchmark = row['benchmark_seq']
    assert len(seq_helix) == len(seq_benchmark)
    identity = 0
    buried_len = 0
    # print(row)
    for i, res in enumerate(seq_helix):
        if i + start in buried_res:
            buried_len += 1
            if res == seq_benchmark[i]:
                identity += 1
    if buried_len == 0:
        print("No buried residues found in chain B")
        print(buried_res)
        return 0

    # print(buried_len)
    # print(100 * identity /buried_len)
    return 100 * (identity / buried_len)


def delta_graph(df, args):
    def relative_to_focusdir(row):
        patchfile = row['patchman_file_x']
        return os.path.join(*patchfile.split('/')[2:])
    df['patchman_basename'] = df.apply(relative_to_focusdir, axis=1)
    prot2 = 'residue_lock_combined_ramp'
    prot1 = 'specialrot_3_combined_ramp'
    df1 = df[df.protocol==prot1]
    df2 = df[df.protocol==prot2]
    dfplot = pd.merge(df1, df2, on=['patchman_basename', 'target',
                                    'name_y'],
                      suffixes=('_{}'.format(prot1),
                                '_{}'.format(prot2)),
                      validate='1:1')
    dfplot = dfplot[dfplot['{}_{}'.format(args['--xaxis'], prot1)] > dfplot['{}_{}'.format(args['--xaxis'], prot2)]]
    dfplot['delta'] = dfplot['{}_{}'.format(args['--yaxis'], prot1)] - dfplot['{}_{}'.format(args['--yaxis'], prot2)]

    bw = 0.2
    xname = '{}_{}'.format(args['--xaxis'], prot1)
    means = []
    hbonds = set(dfplot[xname])
    for n in hbonds:
        mean = dfplot[dfplot[xname]==n]['delta'].mean()
        means.append(mean)
    sns.violinplot(data=dfplot, x=xname, y='delta',
                bw=bw)
    plt.scatter(x=range(len(means)), y=means, c='k', marker='_', s=200)
    plt.tight_layout()
    plt.show()


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    global workspace
    workspace = ws.workspace_from_dir(args['<workspace>'])
    plot_type = args['<plot_type>']
    if args['--test']:
        dfpath = os.path.join(workspace.rifdock_outdir,
                              'combined_benchmark_rev',
                              'test_dataframe.pkl')
    else:
        dfpath = os.path.join(workspace.rifdock_outdir,
            'combined_benchmark_rev', 'final.pkl')
    # dfpath = os.path.join(workspace.root_dir, '..', 'benchmark_analysis', 'final.pkl')
    df_temp = utils.safe_load(dfpath)
    # df_temp = df_temp[df_temp['minimized']]

    if args['--filter'] and not args['--filter-plot']:
        import yaml
        filter_path = args['--filter']
        with open(filter_path, 'r') as f:
            filters = yaml.load(f.read())
        scorelist = []
        for name, group in df_temp.groupby(['patch_len', 'name_x', 'target']):
            scorelist.append(filter.parse_filter(filters, group))
        df_temp = pd.concat(scorelist)

    if args['--filter-plot']:
        df_bench = os.path.join(
            workspace.rifdock_outdir, 'combined_benchmark_rev', 'bench_interfaces_analysis.pkl'
        )
        df_bench = utils.safe_load(df_bench)
        df_bench['benchmark'] = True
        df_temp['benchmark'] = False
        df_temp = pd.concat([df_temp, df_bench])

    if args['--id-cutoff']:
        id_cutoff_path = os.path.join(workspace.rifdock_outdir,
                'combined_benchmark_rev', 'below_{}_pid.pkl'.format(
                    args['--id-cutoff']
                    ))
        if not os.path.exists(id_cutoff_path):
            id_cutoff = float(args['--id-cutoff'])
            print('Dataframe initial size: {}'.format(df_temp.shape[0]))
            df = pd.DataFrame()
            for name, group in df_temp.groupby(['name_x', 'chain']):
                pdbid = name[0]
                chain = name[1]
                homologs = homology.find_homologs(pdbid, id_cutoff,
                        chain=chain)
                group['patchman_pdbid'] = group.apply(get_patchman_pdbid,
                        axis=1)
                group = group[~group['patchman_pdbid'].isin(homologs)]
                df = pd.concat([df, group])
            df.to_pickle(id_cutoff_path)
        else:
            df = utils.safe_load(id_cutoff_path)
        print('Dataframe final size: {}'.format(df.shape[0]))
    else:
        df = df_temp

    if args['--length']:
        length_name = workspace.scaffold_prefix + args['--length']
        df = df[df['patch_len'] == length_name]

    if args['--focus-dir']:
        focusdirs = args['--focus-dir'].split(',')
        df = df[df['focus_dir'].isin(focusdirs)]
    print('Possible values for x- and y-axes are:')
    for col in df.columns:
        print(col)
    if args['--yaxis'] == 'helix_percent_identity' or args['--hue'] ==\
            'helix_percent_identity':
        df['helix_percent_identity'] = df.apply(calc_identity, axis=1)
        df.to_pickle('temp.pkl')
        df = pd.read_pickle('temp.pkl')
    if not float(args['--rmsd-cutoff']) == 0:
        df = df[df['best_rmsd'] < float(args['--rmsd-cutoff'])]

    if args['--yaxis'] == 'hydrophobic_identity':
        df['hydrophobic_identity'] = df.apply(calc_hydrophobic_identity,
                axis=1)

    if args['--yaxis'] == 'buried_identity':
        buried_id_df_name = os.path.join(workspace.rifdock_outdir, 
                'combined_benchmark_rev', 
                'buried_identity_{}A'.format(args['--rmsd-cutoff']) +\
                '.pkl')
        if not os.path.exists(buried_id_df_name):
            df['buried_len'] = df.apply(calc_buried_identity, axis=1)
            df.to_pickle(buried_id_df_name)
        else:
            df = utils.safe_load(buried_id_df_name)

    design_id = False
    if args['--yaxis'] == 'design-identity' or \
            args['--yaxis'] == 'design_identity' or \
            args['--yaxis'] == 'design-id':
        design_id = True
        design_id_type = '--yaxis'
    if args['--xaxis'] == 'design-identity' or \
            args['--xaxis'] == 'design_identity' or \
            args['--xaxis'] == 'design-id':
        design_id = True
        design_id_type = '--xaxis'
    if design_id:
        df[args[design_id_type]] = df.apply(calc_designed_identity,
                axis=1)
        df = df[pd.notna(df[args[design_id_type]])]

    if args['--patch-id-cutoff']:
        print('Filtering on patch ID. Dataframe initial size: {}'.format(df.shape[0]))
        df = df[df['sequence_identity'] < float(args['--patch-id-cutoff'])]
        print('Dataframe final size: {}'.format(df.shape[0]))
    if not args['--pose-score-cutoff'] == 'False':
        df = df[df['pose_score'] < float(args['--pose-score-cutoff'])]

    if args['--target']:
        df = df[(df['name_x'] == args['--target'].split('_')[0]) &
                (df['target'] == args['--target'].split('_')[1]
            )]

    if args['--protocol']:
        df = df[df.protocol==args['--protocol']]

    if plot_type=='logo':
        import weblogo
        df['bench_start'] = df.apply(lambda x: x['rosetta_resis'][0],
                axis=1)
        for name, group in df.groupby([args['--hue'], 'chain',
            'bench_start', 'benchmark_resis']):
            create_web_logo(name, group)

    global y_names
    y_names = {
        'contact_molecular_surface': 'Contact molecular surface',
        'buns_all': 'Interface buried unsatisfied hydrogen bonds',
        'n_hbonds': 'Cross-interface hydrogen bonds',
        'contact_molecular_surface_perres': 'Per-residue contact\nmolecular surface',
        'n_hbonds_perres': 'Per-residue interface\nhydrogen bonds/BUNS',
        'buns_all_perres': 'Per-residue interface buried \nunsatisfied hydrogen bonds',
        'value': 'Per-residue #\nH-bonds or BUNS',
        'interface_dG_perres': 'Per-residue interface ΔG',
    }
    if args['--xaxis'] == 'protocol':
        global order
        global labels
        order = ['base',
                 'buns_penalty', #'buns_penalty_pruned',
                 #'deleteme',
                 'specialres', 'combined',
                 'residue_lock', 'residue_lock_combined', 'residue_lock_combined_ramp',
                 'specialrot',
                 'specialrot_combined', 'special_combined_ramp', #'combined_nomin',
                 'specialrot_3', 'specialrot_3_combined', 'specialrot_3_combined_ramp',] #'interface']
        labels = ['Base',
                  'BUNS penalty', #'BUNS pen. pruned',
                  #'NativeResidue',
                  'Special res.', '+ BUNS',
                  'Res. lock',  '+ BUNS', '+ ramp cst.',
                  'Special rot.',
                  '+ BUNS', '+ ramp cst.', #'Special rot. + BUNS\n(no cst)',
                  'Special rot. (wt. 3)', '+ BUNS', '+ ramp cst.',] #'INTERFACE']
        color_dict = {
            'base': palette['darkgray'],
            'buns_penalty': palette['yellow'],
            'buns_penalty_pruned': palette['yellow'],
            'deleteme': palette['red'],
            'specialres': palette['red'],
            'combined': palette['red'],
            'residue_lock': palette['green'],
            'residue_lock_combined': palette['green'],
            'residue_lock_combined_ramp': palette['green'],
            'specialrot': palette['purple'],
            'specialrot_combined': palette['purple'],
            'special_combined_ramp': palette['purple'],
            'combined_nomin': palette['purple'],
            'specialrot_3': palette['blue'],
            'specialrot_3_combined': palette['blue'],
            'specialrot_3_combined_ramp': palette['blue'],
            'interface' : '#A5AA99',
        }
        global colors
        colors = []
        for prot in order:
            colors.append(color_dict[prot])
            # if args['--buns-vs-bonds']:
            #     colors.append(palette['lightgray'])
        sns.set_palette(sns.color_palette(colors))

    
    if plot_type == 'seq':
        plot_sequence_recovery(df, args)
    if plot_type == 'scatter':
        scatterplot(df, args) 
    if plot_type == 'bar':
        barplot(df, args)
    if plot_type == 'pvp':
        protocol_vs_protocol(df, args)
    if plot_type == 'delta':
        delta_graph(df, args)


if __name__=='__main__':
    main()
