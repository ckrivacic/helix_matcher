'''
Usage:
    prep_rifgen.py <target_pdb> <output_folder> [options]

options:
    --chain=LETTER, -c  
    Which chain of the protein to use
TO DO:
    Allow user to input a chain, figure out what Rosetta chain that is,
    and split that chain
    Allow user to specify ranges
'''
from patches import Patches
import docopt
from pyrosetta import pose_from_file
from pyrosetta import init
import numpy as np
import sys, os


def write_flags(folder, tarpath):
    flags_rifgen = '''
-rifgen:data_cache_dir .

-database {db}

-rifgen::rif_type RotScore  # use RotScoreSat for small molecules

-rifgen:target      {target}
-rifgen:target_res residues.txt
-rifgen:outdir     .
-rifgen:outfile    test_out.rif.gz

-rifgen:score_cut_adjust 0.8

-rifgen::rosetta_field_resl 0.25          # grid spacing for vdw/lksol scoring
-rifgen::search_resolutions 3.0 1.5 0.75  # search resl for denovo hydrophobic interactions

# same name as rif_dock_test flag, but not same!
# this roughly controls how many apolar rotamers are 'docked' for denovo apolar rif gen
# default 1 billion
-rifgen:beam_size_M 1000.0

# dump for sanity
-rifgen:rif_hbond_dump_fraction 0.01
-rifgen:rif_apo_dump_fraction 0.01

-rifgen:extra_rotamers false
-rifgen:extra_rif_rotamers true

# rosetta stuff
-renumber_pdb
-add_orbitals

-hbond_cart_sample_hack_range 0.33
-hbond_cart_sample_hack_resl  0.33
-rifgen:tip_tol_deg        60.0 # for now, do either 60 or 36
-rifgen:rot_samp_resl       6.0

-rifgen:hash_preallocate_mult 0.125
-rifgen:max_rf_bounding_ratio 4.0

# geometry of bounding grids
-rifgen:hash_cart_resls   16.0   8.0   4.0   2.0   1.0
-rifgen:hash_cart_bounds   512   512   512   512   512
-rifgen:lever_bounds      16.0   8.0   4.0   2.0   1.0
-rifgen:hash_ang_resls     38.8  24.4  17.2  13.6  11.8 # yes worky worky
-rifgen:lever_radii        23.6 18.785501 13.324600  8.425850  4.855575

    '''.format(db='/wynton/home/kortemme/krivacic/software/rosetta_rifdock/database', target=tarpath)

    if not os.path.exists(folder):
        os.mkdir(folder)
    with open(os.path.join(folder, 'flags'), 'w') as f:
        f.write(flags_rifgen)



def write_to_file(reslist, folder):
    f = open(os.path.join(folder, 'residues.txt'), 'w')
    for res in reslist:
        f.write(str(res) + '\n')


def test_3n2n():
    init()
    args = docopt.docopt(__doc__)
    # posefile = os.path.abspath('../test_files/3n2n.pdb')
    posefile = os.path.abspath(args['<target_pdb>'])
    pose = pose_from_file(posefile)
    patches = Patches(pose)
    ranges = [(50,138), (150, 173), (201, 220)]
    reslist = []
    for r in ranges:
        reslist.extend(np.arange(r[0]-35, r[1]-35))
    patches.set_reslist(reslist)
    print(patches.reslist)
    patches.determine_surface_residues()
    patches.map_residues()

    parent_folder = os.path.abspath(os.path.join(args['<output_folder>']))
    target_pdb = os.path.join(parent_folder, 'target.pdb')
    i = 1
    for res in patches.reslist:
        print('RES CENTER {}'.format(res))
        patch_folder = os.path.join(parent_folder, 'patch_{}'.format(i))
        i += 1
        if not os.path.exists(patch_folder):
            os.makedirs(patch_folder, exist_ok=True)
        print(patches.nearest_n_residues(res, 100, cutoff=10.5,
            pymol=True))
        write_to_file(patches.nearest_n_residues(res, 100, cutoff=10.5),
                patch_folder)
        write_flags(patch_folder, target_pdb)

    pose.dump_pdb(target_pdb)

    os.symlink(os.environ['RIFGEN'], os.path.join(
        parent_folder, 'rifgen'
        ))
    os.symlink(os.environ['RIFDOCK'], os.path.join(
        parent_folder, 'rifdock'
        ))


def main():
    init()
    args = docopt.docopt(__doc__)
    posefile = os.path.abspath(args['<target_pdb>'])
    pose = pose_from_file(posefile)#.split_by_chain(1)

    if args['--chain']:
        print('MAKING PATCHES FOR CHAIN {}'.format(args['--chain']))
        for i in range(1, pose.num_chains()):
            info = pose.split_by_chain(i).pdb_info().pose2pdb(1)
            if info.split(' ')[1] == args['--chain']:
                pose = pose.split_by_chain(i)
                if pose.size() < 5:
                    raise('Error: chain too small.')
                break
    else:
        pose = pose.split_by_chain(1)
    patches = Patches(pose)
    patches.determine_surface_residues()
    print(patches.reslist)
    patches.map_residues()
    print(patches.resmap)
    parent_folder = os.path.abspath(os.path.join(args['<output_folder>']))
    target_pdb = os.path.join(parent_folder, 'target.pdb')
    i = 1
    for res in patches.reslist:
        patch_folder = os.path.join(parent_folder, 'patch_{}'.format(i))
        i += 1
        if not os.path.exists(patch_folder):
            os.makedirs(patch_folder, exist_ok=True)
        print(patches.nearest_n_residues(res, 100, cutoff=10.5,
            pymol=True))
        write_to_file(patches.nearest_n_residues(res, 100, cutoff=10.5),
                patch_folder)
        write_flags(patch_folder, target_pdb)

    pose.dump_pdb(target_pdb)

    os.symlink(os.environ['RIFGEN'], os.path.join(
        parent_folder, 'rifgen'
        ))
    os.symlink(os.environ['RIFDOCK'], os.path.join(
        parent_folder, 'rifdock'
        ))


if __name__=='__main__':
    # main()
    test_3n2n()
