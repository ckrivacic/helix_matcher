"""
Runs PatchMAN.

Usage:
    patchman.py <workspace> [options]

Options:
    --target=PDB, -t  Only runs on a specific target
    --task=INT  Only run a specific task
    --flexpepdock  Run flexpepdock on matched motifs
    --relax  Do not run relax on target
"""
from pyrosetta import pose_from_pdb
from pyrosetta import init
from pyrosetta import pose_from_file
from distutils.dir_util import copy_tree
import docopt
import os, sys, glob
from helix import workspace as ws
from helix.utils import utils
from helix import big_jobs
from helix.patchman import extract_peps_for_motif
from Bio.SubsMat import MatrixInfo


def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector


def score_sequences(seq1, seq2):
    matrix = MatrixInfo.blosum62
    score = 0
    assert(len(seq1)==len(seq2))
    for i in range(0, len(seq1)):
        aa1 = seq1[i]
        aa2 = seq2[i]
        if (aa1, aa2) in matrix:
            score += matrix[(aa1, aa2)]
        else:
            score += matrix[(aa2, aa1)]

    return score



def align_matches(folder, matches, workspace, patch):
    '''
    For each match, get the sequence of the patch and the match. Record
    sequence compatibility.
    '''
    import ast
    import pandas as pd
    # from Bio import pairwise2
    # from Bio.pairwise2 import format_alignment

    # init('-ignore_zero_occupancy false')
    init()
    dict_list = []
    print('Evaluating sequence similarity for {}'.format(
        matches 
        ))
    split = patch.split('_')
    patchno = split[0]
    basename = split[1]
    chain = split[2]
    patch_pdb = "{}_{}.pdb".format(
            patchno, basename, chain
            )
    patch_pose = pose_from_file(patch_pdb)
    patch_sequence = patch_pose.sequence()

    print('FILES IN FOLDER')
    print(os.listdir())

    latest_pdbid = None
    nline = 0
    line_idx = 0
    for match in matches:
        nline += 1
        position_list = match[7]
        # filename = line.strip().split(' ')[1].split('/')[-1]
        # match_pdbid = filename.split('_')[0]
        # match_chain = filename.split('.')[0].split('_')[1]
        match_pdbid = match[1]
        match_chain = match[6]

        complexes = []
        # if match_pdbid.lower() != '1m6y':
            # continue
        # print("GLOBSTR")
        # print("{}_{}_{}_*.pdb".format(patchno,
            # match_pdbid.lower(), line_idx))
        for comp in glob.glob('{}_{}_{}_*.pdb'.format(patchno,
            match_pdbid.lower(), line_idx)):
            complexes.append(os.path.join(
                folder, 'docked_full', comp
                ))
        if len(complexes) < 1:
            # print('No complexes for this match; skipping')
            line_idx += 1
            continue
        else:
            print('PDBID')
            print(match_pdbid)
            print('The following complexes were found:')
            print(complexes)

        # All this try/except nonsense is probably not necessary
        if not match_pdbid == latest_pdbid:
            try:
                match_pose = utils.pose_from_wynton(match_pdbid)
            except:
                try:
                    if match_pdbid.lower() == '4k0f':
                        match_pose = utils.pose_from_wynton('5eqb')
                except:
                    print('Could not find PDB {}'.format(match_pdbid))
                    continue
        match_pose = utils.pose_get_chain(match_pose, match_chain)
        match_sequence = ''
        for pos in position_list:
            print(pos)
            match_sequence += match_pose.sequence(pos[0]+1,
                    pos[1]+1) 
        # MASTER output is apparently 0-indexed

        print('Aligning the following sequences')
        print(patch_sequence)
        print(match_sequence)
        score = score_sequences(patch_sequence, match_sequence)
        print('SCORE: {}'.format(score))
        # alignments = pairwise2.align.globalds(patch_sequence,
                # match_sequence, matrix, -1, -0.3)
        # print(format_alignment(*alignments[0]))
        # print(alignments[0].score)
        for cmplx in complexes:
            pdb_save = os.path.relpath(cmplx,
                    start=workspace.root_dir)
            dict_list.append({
                'patch': patchno,
                'complex': pdb_save,
                'target_pdb': basename,
                'patch_pdb': os.path.join(folder, patch_pdb),
                'patch_sequence': patch_sequence,
                'match_pdb': match_pdbid,
                'match_chain': match_chain,
                'match_resis': position_list,
                'match_sequence': match_sequence,
                'alignment_score': score,
                })
        print('Finished {} lines'.format(nline))
        line_idx += 1
        # if match_pdbid.lower() == '4m8r':
        # if match_pdbid.lower() == '1m6y':
            # print(dict_list)
            # sys.exit()

    print('Got through {} lines total'.format(line_idx))

    return pd.DataFrame(dict_list)


def main():
    args = docopt.docopt(__doc__)
    workspace, job_info = big_jobs.initiate()
    if not hasattr(workspace, 'docking_directory'):
        raise Exception("Error: run_sge.py requires RIFWorkspaces as input. You "\
                "may have provided the root directory to this script "\
                "somehow.")

    total_jobs = len(glob.glob(workspace.focus_dir + '/patch_*/len_*'))
    print('TOTAL JOBS: {}'.format(total_jobs))
    num_tasks = total_jobs

    task_id = job_info['task_id']
    if not task_id:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0

    print('TASK: {}'.format(task_id))

    if 'TMPDIR' in os.environ:
        os_tmp = os.environ['TMPDIR']
    else:
        os_tmp = os.path.join(os.environ['HOME'], 'temp')
    outdir_temp = os.path.join(os_tmp, str(task_id))
    if not os.path.exists(outdir_temp):
        os.makedirs(outdir_temp, exist_ok=True)

    folders = workspace.patches
    folder = folders[task_id]

    print("Running PatchMAN for {}".format(folder))
    basename = os.path.basename(folder)
    tempdir = os.path.join(outdir_temp, basename)
    print('TEMPDIR IS {}'.format(tempdir))
    copy_tree(folder, tempdir)
    os.chdir(tempdir)

    motif_path = os.path.join(tempdir, '???_target.pdb')
    motif_list = glob.glob(motif_path)
    if len(motif_list) > 1:
        print('ERROR: More than one input motif found. Motif list:')
        print(motif_list)
        sys.exit()
    else:
        motif = motif_list[0]

    target = workspace.target_path_clean
    print("TARGET:")
    print(target)

    # Split surface into structural patches (done in prep_patchman)
    # print('Running split_to_motifs.py on {}'.format(target))
    # cmd = [os.path.join(workspace.patchman_path, 'split_to_motifs.py'),
            # target]
    # utils.run_command(cmd)

    # Create a list of motifs
    os.system("ls ???_target.pdb > motif_list")

    # Create pds files for running MASTER
    cmd = [os.path.join(workspace.master_path, 'createPDS'), '--type',
            'query', '--pdbList', 'motif_list']
    utils.run_command(cmd)

    # Create a list of database structures for the template search
    os.system("ls {} > db_list".format(
        os.path.join(workspace.master_db_path, '*', '*pds')
        ))

    # Remove homologs
    homologs = []
    if os.path.exists(os.path.join(workspace.focus_dir, 'homologs')):
        with open(os.path.join(workspace.focus_dir, 'homologs'), 'r') as f:
            for line in f:
                homologs.append(line.strip().lower())
    matches = glob.glob('*_matches')
    fout = open('db_list_nohom', 'w')
    removed = open('removed_psds', 'w')
    with open('db_list', 'r') as f:
        for line in f:
            pdbid = line.split('/')[-1].split('_')[0]
            if pdbid.lower() != os.path.basename(target).split('_')[0].lower():
                fout.write(line)
            else:
                removed.write(line)
    fout.close()
    removed.close()
    os.remove('db_list')
    os.rename('db_list_nohom', 'db_list')

    # Run MASTER for all motifs
    pds_list = glob.glob('*pds')
    for pds in pds_list:
        name = pds.split('.')[0]
        exe = os.path.join(workspace.master_path, 'master')
        cmd = [exe, '--query', pds,
                '--targetList', 'db_list',
                '--bbRMSD', '--rmsdCut', '1.5', '--topN', '1000000',
                '--matchOut', '{}_matches'.format(name)]
        utils.run_command(cmd)

    # Prepack the receptor structure for further FlexPepDock refinement
    if args['--relax']:
        import platform
        ostype = platform.system()
        if ostype == 'Linux':
            suffix = 'linuxgccrelease'
        elif ostype == 'Darwin':
            suffix = 'macosclangrelease'
        exe = os.path.join(
                workspace.rosetta_dir, 'source', 'bin',
                'relax.{}'.format(suffix)
                )
        cmd = [exe, '-s', target, '-out:pdb', '-scorefile', 'ppk.score.sc',
                '-nstruct', '1', '-ex1', '-ex2aro', '-use_input_sc']
        utils.run_command(cmd)

        os.system("grep ' A ' *target.clean_0001.pdb > target.ppk.pdb")
    else:
        os.system("grep ' A ' {} > target.ppk.pdb".format(target))

    # Extract templates and thread the peptide sequence
    length = os.path.basename(folder).split('_')[-1]
    all_matches = []
    init()
    with open('motif_list', 'r') as f:
        for line in f:
            line = line.strip('\n')
            name = line.split('.')[0]
            # script = os.path.join(
                    # workspace.patchman_path,
                    # 'extract_peps_for_motif.py'
                    # )
            # cmd = [workspace.python_path, script, '-m', name + '_matches', '-d', '-r',
                    # 'target.ppk.pdb', '--patch', line, '-l', str(length)]
            # utils.run_command(cmd)
            arglist = ['-m', name + '_matches', '-d', '-r',
                    'target.ppk.pdb', '--patch', line, '-l',
                    str(length)]
            motif_args = extract_peps_for_motif.arg_parser().parse_args(arglist)
            # motif_args.match_list = name + '_matches'
            # motif_args.design = True
            # motif_args.receptor = 'target.ppk.pdb'
            # motif_args.patch = line
            # motif_args.peplen = str(length)
            # print('ARGLIST')
            # print(arglist)
            # motif_args = ' '.join(arglist)
            print('MOTIF ARGS')
            print(motif_args)
            matches = extract_peps_for_motif.main(motif_args)

            # Create a list of input structures for refinement

            alignment_df = align_matches(tempdir, matches,
                    workspace, line)

    os.makedirs('docked_full/', exist_ok=True)
    alignment_df.to_pickle('alignment_scores.pkl')

    # Run FlexPepDock refinement
    if args['--flexpepdock']:
        from pyrosetta.rosetta.core.pack.task import TaskFactory
        from pyrosetta.rosetta.core.pack.task import operation
        from pyrosetta.rosetta.core.select import residue_selector
        from pyrosetta import create_score_function
        import pyrosetta

        os.system('ls ???_????_*_*.pdb > input_list')
        inputs = []
        with open('input_list', 'r') as f:
            for line in f:
                inputs.append(line.strip())
        for pdb in inputs:
            pose = pose_from_pdb(pdb)
            ref = create_score_function('ref2015')
            fastdes = pyrosetta.rosetta.protocols.denovo_designs.movers.FastDesign(ref)

            selector = residue_selector.ChainSelector('B')
            not_selector = residue_selector.NotResidueSelector(selector)
            tf = TaskFactory()
            no_packing = operation.PreventRepackingRLT()
            static = operation.OperateOnResidueSubset(no_packing,
                    not_selector)
            notaa = operation.ProhibitSpecifiedBaseResidueTypes(
                    strlist_to_vector1_str(['GLY']),
                    selector)
            tf.push_back(static)
            tf.push_back(notaa)
            packertask = tf.create_task_and_apply_taskoperations(pose)
            print('REPACK')
            print(packertask.repacking_residues())
            print('DESIGN')
            print(packertask.designing_residues())

            fastdes.set_task_factory(tf)
            fastdes.set_scorefxn(ref)
            fastdes.apply(pose)
            score = ref(pose)
            pose.dump_pdb(pdb)

        exe = os.path.join(
                workspace.rosetta_dir, 'source', 'bin',
                'FlexPepDocking.{}'.format(suffix)
                )
        if not os.path.exists('docked_full/'):
            os.path.makedirs('docked_full', exist_ok=True)
        cmd = [exe, '-in:file:l', 'input_list', '-scorefile',
                'docked_full/score.sc',
                '-out:pdb_gz', '-lowres_preoptimize',
                '-flexPepDocking:pep_refine',
                '-out:prefix', 'docked_full/',
                '-flexPepDocking:flexpep_score_only', '-ex1', '-ex2aro',
                '-use_input_sc', 'unboundrot', target]
        utils.run_command(cmd)
        os.system('mv alignment_scores.pkl docked_full/')
    else:
        os.system('mv alignment_scores.pkl docked_full/')
        os.system('mv ???_*_*_*.pdb docked_full/')
        os.system('mv db_list docked_full/')
        os.system('mv removed_psds docked_full/')
        os.system('mv *_matches docked_full/')
        os.system('gzip docked_full/*.pdb')

    # Copy back to main folder
    # os.system('tar -cf docked_full.tar docked_full/')
    outputs = os.path.join(tempdir, 'docked_full/')
    final_out = os.path.join(folder, 'docked_full/')
    copy_tree(outputs, final_out)

if __name__=='__main__':
    main()
