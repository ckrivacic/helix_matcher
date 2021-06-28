'''
Run rifdock, either locally or on SGE.

Usage:
    run_rifdock.py <folder> [Options]

Options:
    --job-distributor=STR, -d  Running on SGE or locally? Defaultse to
    local.
    --tasks=NUM, -n  How many jobs to split into if running on SGE.
    --sge  Run on SGE job distributor
'''
import os, sys, docopt
import glob
import math
import subprocess
from subprocess import Popen, PIPE


def write_flags(folder):
#-rif_dock:target_res            residue_numbers.txt
    
    scaffold = '/home/ckrivacic/software/helix_matcher/test_files/scaffolds.txt'
    tarpath, cache = get_flag_params(folder)
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
    '''.format(db=os.path.join(os.environ['ROSETTA'], 'database'),
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



def main():
    args = docopt.docopt(__doc__)
    print(args)
    folder = os.path.abspath(args['<folder>'])

    total_jobs = len(glob.glob(folder + '/patch_*'))
    if '--tasks' in args:
    #if args['--tasks']:
        num_tasks = int(args['--tasks'])
    else:
        num_tasks = total_jobs

    if '--sge' in args:
    #if args['--sge']:
        task = int(os.environ['SGE_TASK_ID']) - 1
    else:
        task = 0
    
    start_job = task * math.ceil((total_jobs / num_tasks))
    stop_job = start_job + math.ceil(total_jobs / num_tasks) - 1

    folders = sorted(glob.glob(folder + '/patch_*'))

    # rifdock = os.environ['RIFDOCK']
    rifdock = os.path.join(folder, 'rifdock')
    for fold in folders[start_job:stop_job+1]:
        print(fold)
        os.chdir(fold)
        write_flags(fold)
        flags = os.path.join(fold, 'dock_flags')
        myenv['LD_LIBRARY_PATH'] = '/wynton/home/kortemme/krivacic/software/anaconda3/lib/'
        process = Popen([rifdock, '@', flags], stdout=PIPE, stderr=PIPE,
                env=myenv)
        # process = Popen([rifdock, '@', flags], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        exit_code = process.wait()
        if exit_code != 0:
            print(stdout)
            print(stderr)
        #output = subprocess.check_output([rifdock, '@', flags])
        #print(output)


if __name__=='__main__':
    main()
    #get_flag_params('../test_files/6r9d/')
