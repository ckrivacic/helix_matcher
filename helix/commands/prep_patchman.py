"""
Preps a RIFWorkspace for PatchMAN.

Usage:
    helix prep_patchman <workspace> [options]

Options:
    --chain=LETTERS, -c
        Which chain(s) of the protein to use
    --target=PDB, -t
        Which target to do the prep for if multiple exist and you only
        want to focus on one.
    --chainmap=YAML, -m
        YAML file that tells this script which chains go with which
        target.

"""

import docopt
import helix.workspace as ws
from helix.utils import utils
from pyrosetta import init
from pyrosetta import pose_from_file
import os, shutil
import yaml
import glob


def main():
    args = docopt.docopt(__doc__)
    init()
    root_workspace = ws.workspace_from_dir(args['<workspace>'])

    if args['--target']:
        targets = [args['--target']]
    else:
        targets = root_workspace.targets
    
    chainmap = None
    if args['--chainmap']:
        with open(args['--chainmap']) as file:
            chainmap = yaml.load(file)

    sizes = [10, 14, 21, 28]

    orig_dir = os.path.abspath(os.getcwd())
    for target in targets:
        os.chdir(orig_dir)
        try:
            workspace = ws.RIFWorkspace(args['<workspace>'], target)
            if os.path.exists(workspace.target_path_clean):
                continue
            workspace.make_dirs()
            pose = pose_from_file(workspace.initial_target_path)
            chain = None
            if chainmap:
                chain = chainmap[workspace.focus_name]
            elif args['--chain']:
                chain = args['--chain']

            if chain:
                if len(chain) > 1:
                    print("Error: cannot specify more than 1 chain for "\
                            " PatchMAN runs.")
                    sys.exit()
                print('MAKING PATCHES FOR CHAIN {}'.format(chain))
                poses = []
                for i in range(1, pose.num_chains() + 1):
                    try:
                        chainpose = pose.split_by_chain(i)
                    except:
                        continue
                    info = chainpose.pdb_info().pose2pdb(1)
                    if info.split(' ')[1] in chain and chainpose.residue(1).is_protein():
                        if chainpose.size() < 5:
                            raise('Error: chain {} too small.'.format(chain))
                        else:
                            poses.append(chainpose)
                pose = poses[0]
                if len(poses) > 1:
                    for chainpose in poses[1:]:
                        pose.append_pose_to_pose(pose, chainpose)

            else:
                pose = pose.split_by_chain(1)

            target_pdb = os.path.abspath(workspace.target_path)
            pose.pdb_info().set_chains('A')
            pose.dump_pdb(target_pdb)

            script_path = os.path.join(workspace.patchman_path,
                    'split_to_motifs.py')
            os.chdir(workspace.focus_dir)
            utils.run_command([workspace.python_path,
                script_path, target_pdb])
            outputs = glob.glob('???_target.pdb')
            for output in outputs:
                patchno = output[:3]
                patch_folder = os.path.join(workspace.focus_dir,
                        'patch_{}'.format(patchno))
                if not os.path.exists(patch_folder):
                    os.makedirs(patch_folder, exist_ok=True)
                for out_length in workspace.matchlen:
                    len_folder = 'len_{}'.format(out_length)
                    os.makedirs(os.path.join(patch_folder, len_folder),
                            exist_ok=True)
                    shutil.copyfile(output, os.path.join(patch_folder,
                        len_folder, os.path.basename(output)))
                os.remove(output)

        except Exception as e:
            import traceback
            print("Error finding patches for {}. Error was:".format(target))
            print(traceback.format_exc())
            print(e)
