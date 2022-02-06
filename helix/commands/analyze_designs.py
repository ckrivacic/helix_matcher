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
'''
import docopt
import os
import pandas as pd
from klab import scripting
from helix.utils import utils
from helix.utils import homology
import helix.workspace as ws
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from pyrosetta import *
from pyrosetta.rosetta.core.select import residue_selector


def first3_res_correct(pose, res, seq):
    first = pose.residue(res).name1() == seq[0]
    second = pose.residue(res+1).name1() == seq[1]
    third = pose.residue(res+2).name1() == seq[2]
    return first and second and third


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
        print('Residues witheld for {}'.format(row['patchman_file']))
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
    sns.stripplot(data=df, x='protocol', y=args['--yaxis'],
            order=order, color='.5', alpha=0.5)
    sns.violinplot(data=df, x='protocol', y=args['--yaxis'],
            order=order)
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
        sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=hue, size=size, sizes=(50,300),)
    else:
        sns.scatterplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
                hue=hue,)
    plt.show()


def barplot(df, args):
    # order = ['1b33_K', '1b33_K_buns_noprune', '1b33_K_buns_penalty',
            # '1b33_K_nativelike', '1b33_K_specialrot']
    # labels = ['Base', 'BUNS penalty', 'BUNS penalty + prune', 
            # 'Freeze good rotamers', 'Special rotamer bonus']
    order=None
    if args['--hue']:
        hue = args['--hue']
    else:
        hue = None
    fig, ax = plt.subplots()
    # sns.stripplot(data=df, x='focus_dir', y=args['--yaxis'],
            # order=order, color='.5', alpha=0.5, ax=ax)
    sns.barplot(data=df, x=args['--xaxis'], y=args['--yaxis'],
            hue=hue, order=order, ax=ax)
    if args['--xaxis'] == 'focus_dir':
        ax.set_xticklabels(labels)
    plt.show()


def get_patchman_pdbid(row):
    return os.path.basename(row['patchman_file']).split('_')[1].upper()

def calc_buried_identity(row):

    init()
    pose = pose_from_file(os.path.join(
        workspace.root_dir, row['patchman_file']))
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
    for i, res in enumerate(seq_helix):
        if i + start in buried_res:
            buried_len += 1
            if res == seq_benchmark[i]:
                identity += 1
    if buried_len == 0:
        print("No buried residues found in chain B")
        print(buried_res)
        return 0

    # return 100 * (identity / buried_len)
    return buried_len


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
        print('Dataframe final size: {}'.format(df.shape[0]))
    else:
        df = df_temp

    if args['--focus-dir']:
        focusdirs = args['--focus-dir'].split(',')
        df = df[df['focus_dir'].isin(focusdirs)]
    print('Possible values for x- and y-axes are:')
    for col in df.columns:
        print(col)
    if args['--yaxis'] == 'helix_percent_identity':
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
            df['buried_identity'] = df.apply(calc_buried_identity, axis=1)
            df.to_pickle(buried_id_df_name)
        else:
            df = utils.safe_load(buried_id_df_name)

    if args['--patch-id-cutoff']:
        print('Filtering on patch ID. Dataframe initial size: {}'.format(df.shape[0]))
        df = df[df['sequence_identity'] < float(args['--patch-id-cutoff'])]
        print('Dataframe final size: {}'.format(df.shape[0]))
    if not args['--pose-score-cutoff'] == 'False':
        df = df[df['pose_score'] < float(args['--pose-score-cutoff'])]

    df['bench_start'] = df.apply(lambda x: x['rosetta_resis'][0],
            axis=1)

    if args['--target']:
        df = df[(df['name_x'] == args['--target'].split('_')[0]) &
                (df['target'] == args['--target'].split('_')[1]
            )]

    if plot_type=='logo':
        import weblogo
        for name, group in df.groupby([args['--hue'], 'chain',
            'bench_start', 'benchmark_resis']):
            create_web_logo(name, group)

    
    if plot_type == 'seq':
        plot_sequence_recovery(df, args)
    if plot_type == 'scatter':
        scatterplot(df, args) 
    if plot_type == 'bar':
        barplot(df, args)


if __name__=='__main__':
    main()
