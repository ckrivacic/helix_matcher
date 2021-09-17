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
from distutils.dir_util import copy_tree
import docopt
import os, sys, glob
from helix import workspace as ws
from helix.utils import utils
from helix import big_jobs


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
    with open(os.path.join(workspace.focus_dir, 'homologs'), 'r') as f:
        for line in f:
            homologs.append(line.strip().lower())
    matches = glob.glob('*_matches')
    fout = open('db_list_nohom', 'w')
    removed = open('removed_psds', 'w')
    with open('db_list', 'r') as f:
        for line in f:
            pdbid = line.split('/')[-1].split('_')[0]
            if pdbid not in homologs:
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
    with open('motif_list', 'r') as f:
        for line in f:
            line = line.strip('\n')
            name = line.split('.')[0]
            script = os.path.join(
                    workspace.patchman_path,
                    'extract_peps_for_motif.py'
                    )
            cmd = [workspace.python_path, script, '-m', name + '_matches', '-d', '-r',
                    'target.ppk.pdb', '--patch', line, '-l', str(length)]
            utils.run_command(cmd)

    # Create a list of input structures for refinement
    os.makedirs('docked_full/', exist_ok=True)

    # Run FlexPepDock refinement
    if args['--flexpepdock']:
        os.system('ls ???_????_*_*.pdb > input_list')
        exe = os.path.join(
                workspace.rosetta_dir, 'source', 'bin',
                'FlexPepDocking.{}'.format(suffix)
                )
        if not os.path.exists('docked_full/'):
            os.path.makedirs('docked_full', exist_ok=True)
        cmd = [exe, '-in:file:l', 'input_list', '-scorefile', 'score.sc',
                '-out:pdb', '-lowres_preoptimize',
                '-flexPepDocking:pep_refine',
                '-out:prefix', 'docked_full/',
                '-flexPepDocking:flexpep_score_only', '-ex1', '-ex2aro',
                '-use_input_sc', 'unboundrot', target]
        utils.run_command(cmd)
    else:
        os.system('mv ???_????_*_*.pdb docked_full/')
        os.system('mv db_list docked_full/')
        os.system('mv removed_psds docked_full/')

    # Copy back to main folder
    os.system('gzip docked_full/*.pdb')
    os.system('tar -cf docked_full.tar docked_full/')
    outputs = os.path.join(tempdir, 'docked_full.tar')
    final_out = os.path.join(folder, 'docked_full.tar')
    copy_tree(outputs, final_out)

if __name__=='__main__':
    main()
