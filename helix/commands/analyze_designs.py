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
        [default: 1.0]

    --pose-score-cutoff=FLOAT, -p  Elimiate poses with a total score
        greater than this float  [default: 0]

    --focus-dir=STR, -f  Only plot for a specific focus dir or list of focus
        dirs (ex. 1b33_K,1b33_K_specialrot)

    --size=COLUMN, -s  Size scatter plot points by this column

    --hue=COLUMN, -h  Color by a column

    --id-cutoff=FLOAT  Filter out results that have a sequence identity
        within this #

    --patch-id-cutoff=FLOAT  filter out patches that have htis amount of
        sequence identity

    --buried-identity  Calculate buried sequence identity. WARNING: Do
        not use without a cutoff... slow.

    --target=STR  Only show results for a specific target

    --protocol=STR  Only show results for this protocol
'''
import docopt
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
    order = None
    order = ['base', #'deleteme', 
            'buns_penalty', 'buns_penalty_pruned',
            'residue_lock', 'specialrot', 'combined',
            'combined_nomin']
    labels = ['Base', #'NativeResidue', 
            'BUNS penalty', 'BUNS pen. pruned', 
            'Residue lock', 'Special rotamer*', 'Combined', 
            'Combined (no cst)']
    # sns.stripplot(data=df, x='protocol', y=args['--yaxis'],
            # order=order, color='.5', alpha=0.5)
    df = df[df['protocol'] != 'deleteme']
    means = df.groupby('protocol')[args['--yaxis']].mean()
    fig, ax = plt.subplots()
    sns.violinplot(data=df, x='protocol', y=args['--yaxis'],
            order=order, )
    plt.scatter(x=range(len(means)), y=means, c='k', marker='_', s=200)
    if args['--xaxis'] == 'protocol':
        ax.set_xticklabels(labels)
        plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


def scatterplot(df, args):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    order = None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    if args['--size']:
        size = args['--size']
        ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=hue, size=size, sizes=(50,300), picker=True)
    else:
        ax = sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=hue, picker=True)
    click = plotting.ClickablePlot(ax, df, args, workspace)

    plt.show()


def barplot(df, args):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    # labels = ['Base', 'BUNS penalty', 'BUNS penalty + prune', 
            # 'Freeze good rotamers', 'Special rotamer bonus']
    if args['--xaxis'] == 'protocol':
        order = ['base', #'deleteme', 
                'buns_penalty', 'buns_penalty_pruned',
                'residue_lock', 'specialrot', 'combined',
                'combined_nomin']
        labels = ['Base', #'NativeResidue', 
                'BUNS penalty', 'BUNS pen. pruned', 
                'Residue lock', 'Special rotamer*', 'Combined', 
                'Combined (no cst)']
    else:
        order=None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    fig, ax = plt.subplots()
    # sns.stripplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
            # order=order, color='.5', alpha=0.5, ax=ax)
    sns.barplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
            hue=hue, order=order, ax=ax)
    if args['--xaxis'] == 'protocol':
        ax.set_xticklabels(labels)
        plt.xticks(rotation=45)
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


def main():
    mpl.use('tkagg')
    args = docopt.docopt(__doc__)
    global workspace
    workspace = ws.workspace_from_dir(args['<workspace>'])
    plot_type = args['<plot_type>']
    dfpath = os.path.join(workspace.rifdock_outdir,
            'combined_benchmark_rev', 'final.pkl')
    df_temp = utils.safe_load(dfpath)

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

    
    if plot_type == 'seq':
        plot_sequence_recovery(df, args)
    if plot_type == 'scatter':
        scatterplot(df, args) 
    if plot_type == 'bar':
        barplot(df, args)
    if plot_type == 'pvp':
        protocol_vs_protocol(df, args)


if __name__=='__main__':
    main()
