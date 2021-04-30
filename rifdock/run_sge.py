'''
Run rifgen on the cluster. One task per patch.
Uses node scratch space to store RIF files created during RIFGen, then
runs RIFDock before copying them back over to the user's patch folder.

Usage:
    cluster_run.py <folder> <scaffold> [options]

Options:
    --sge  Running on the cluster?
    --task=NUM  Which task number, for testing
'''

import sys, os, docopt
import glob
import math
from subprocess import Popen, PIPE, STDOUT
from io import StringIO
from shutil import copytree
from shutil import copyfile
from distutils.dir_util import copy_tree



def write_flags(folder, scaffold):
#-rif_dock:target_res            residue_numbers.txt
    
    tarpath, cache = get_flag_params(folder)
    scaffold = os.path.abspath(scaffold)
    rosetta_path = '/wynton/home/kortemme/krivacic/software/rosetta_rifdock/'
    flags_rifgen = '''
-rif_dock:target_pdb            ./{target}.rif.gz_target.pdb.gz
-rif_dock:target_rf_resl        0.25
-rif_dock:target_rf_cache       ./{cache}
-rif_dock:target_bounding_xmaps ./{target}.rif.gz_BOUNDING_RIF_16.xmap.gz
-rif_dock:target_bounding_xmaps ./{target}.rif.gz_BOUNDING_RIF_08.xmap.gz
-rif_dock:target_bounding_xmaps ./{target}.rif.gz_BOUNDING_RIF_04.xmap.gz
-rif_dock:target_bounding_xmaps ./{target}.rif.gz_BOUNDING_RIF_02.xmap.gz
-rif_dock:target_bounding_xmaps ./{target}.rif.gz_BOUNDING_RIF_01.xmap.gz
-rif_dock:target_rif            ./{target}.rif.gz
-rif_dock:extra_rotamers        0
-rif_dock:extra_rif_rotamers    1
-rif_dock:rot_spec_fname        ./rotamer_index_spec.txt
-database {db}
-rif_dock:rotrf_cache_dir cache/
-rif_dock:data_cache_dir  .
-rif_dock:cache_scaffold_data true  # set to false if you don't want to use/generate

-rif_dock:scaffolds {scaffold} # list of scaffold pdb files

 # optional list of 'res' files for each scaffold
 # must be either one file which applies to all input scaffolds
 # or a list exactly matching the number of input scaffolds
-rif_dock:scaffold_res


# this is where the output will go, and how much
-rif_dock:outdir docked_full/
-rif_dock:dokfile all.dok
-rif_dock:n_pdb_out 20 # max number of output pdbs

# optional flag to add extra output file tag
#-rif_dock:target_tag conf01


# set to true to align all output to scaffolds instead of target
# mostly useful for small molecule binder design
-rif_dock:align_output_to_scaffold false

# include some pikaa lines in output pdbs for rif residues(??)
-rif_dock:pdb_info_pikaa false # this is default I think

# *IFF* using a rif with satisfaction constraints (e.g. RotScoreSat)
# set this to the number of hbonds to the target which are required
# no way yet to explicity say which hbonds are required, but this gives
# some control. searches will be much faster if this number is higher
# of course, don't make it higher than the number of hbonds that can
# actually be made
-require_satisfaction 0
    '''.format(db=os.path.join(rosetta_path, 'database'),
            target=tarpath, cache=cache, scaffold=scaffold)

    if not os.path.exists(folder):
        os.mkdir(folder)
    with open(os.path.join(folder, 'dock_flags'), 'w') as f:
        f.write(flags_rifgen)


def get_flag_params(folder):
    import re
    cache_pattern = '__RF_\w+.pdb_CEN_trhash\d+_resl\d+\.\d+_osamp2_replonlybdry'
    target_pattern = '\w+.rif.gz_target.pdb.gz'
    for f in glob.glob(folder + '/*.gz'):
        cachematch = re.findall(cache_pattern, f)
        targetmatch = re.findall(target_pattern, f)
        if len(cachematch) > 0:
            rfcache = cachematch[0]
        if len(targetmatch) > 0:
            target = targetmatch[0]

    target = os.path.basename(target).split('.')[0]
    return target, rfcache


def execute(cmd, environment=None):
    print('Running command:')
    print(' '.join(cmd))
    # if args['--task']:
    # task = int(args['--task'])
    # else:
    jobno = os.environ['JOB_ID']
    task = os.environ['SGE_TASK_ID']
    logdir = '/wynton/home/kortemme/krivacic/intelligent_design/helix_matcher/rifdock/logs/'
    logfile_path = os.path.join(logdir, 'rifdock_full.o{}.{}'.format(jobno,
        task))

    logfile = open(logfile_path, 'a')
    if not environment:
        environment = os.environ.copy()
    import time

    with Popen(cmd, stdout=PIPE, stderr=STDOUT, bufsize=1,
            env=environment,
            universal_newlines=True) as p, StringIO() as buf:
        for line in p.stdout:
            print(line, end='')
            buf.write(line)
            # logfile.write(line)
            # sys.stdout.buffer.write(line)
        output = buf.getvalue()
    rc = p.returncode

    return rc


def run_command(cmd, environment=None):
    print("Working directory: {}".format(os.getcwd()))
    print("Command: {}".format(' '.join(cmd)))
    sys.stdout.flush()

    process = Popen(cmd, env=environment)

    print("Process ID: {}".format(process.pid))
    print()
    sys.stdout.flush()

    process.wait()


def main():
    args = docopt.docopt(__doc__)
    folder = os.path.abspath(args['<folder>'])

    total_jobs = len(glob.glob(folder + '/patch_*'))
    print('TOTAL JOBS: {}'.format(total_jobs))
    num_tasks = total_jobs

    if args['--sge']:
        task = int(os.environ['SGE_TASK_ID']) - 1
    else:
        if args['--task']:
            task = int(args['--task'])
        else:
            task = 0

    print('TASK: {}'.format(task))

    start_job = task * math.ceil(total_jobs / num_tasks)
    stop_job = start_job + math.ceil(total_jobs / num_tasks)
    print('START JOB: {}'.format(start_job))
    print('STOP JOB: {}'.format(stop_job))

    folders = sorted(glob.glob(folder + '/patch_*'))

    rifgen = os.path.join(folder, 'rifgen')
    rifdock = os.path.join(folder, 'rifdock')

    if 'TMPDIR' in os.environ:
        os_tmp = os.environ['TMPDIR']
    else:
        os_tmp = os.path.join(os.environ['HOME'], 'temp')
    outdir_temp = os.path.join(os_tmp, str(task))
    if not os.path.exists(outdir_temp):
        os.makedirs(outdir_temp, exist_ok=True)

    target = os.path.join(folder, 'target.pdb')
    new_target = os.path.join(outdir_temp, 'target.pdb')
    copyfile(target, new_target)

    for fold in folders[start_job:stop_job]:
        print('Running RIFGEN for {}'.format(fold))
        basename = os.path.basename(fold)
        tempdir = os.path.join(outdir_temp, basename)
        print('TEMPDIR IS {}'.format(tempdir))
        copy_tree(fold, tempdir)
        os.chdir(tempdir)
        flags = os.path.join(tempdir, 'flags')
        myenv = os.environ.copy()
        myenv['LD_LIBRARY_PATH'] = '/wynton/home/kortemme/krivacic/software/anaconda3/lib/'

        # exit_code = execute([rifgen, '@', flags], environment=myenv)
        run_command([rifgen, '@', flags], environment=myenv)
        # if exit_code != 0:
            # print(stdout)
            # print(stderr)

        # write_flags(fold)
        print('Prepping RIFDOCK for {}'.format(fold))
        write_flags(tempdir, args['<scaffold>'])

        flags = os.path.join(tempdir, 'dock_flags')
        print('Running RIFDOCK for {}'.format(fold))
        # exit_code = execute([rifdock, '@', flags], environment=myenv)
        run_command([rifdock, '@', flags], environment=myenv)
        # if exit_code != 0:
            # print(stdout)
            # print(stderr)

        # else:
        outputs = os.path.join(tempdir, 'docked_full')
        final_out = os.path.join(fold, 'docked_full')
        copy_tree(outputs, final_out)
        

if __name__ == '__main__':
    main()
