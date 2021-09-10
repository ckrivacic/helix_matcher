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


def main():
    args = docopt.docopt(__doc__)
    root_workspace = ws.workspace_from_dir(args['<workspace>'])

    if args['--target']:
        targets = args['--target']
    else:
        targets = root_workspace.targets
    
    chainmap = None
    if args['--chainmap']:
        with open(args['--chainmap']) as file:
            chainmap = yaml.load(file)

    sizes = [10, 14, 21, 28]

    for target in targets:
        try:
            workspace = ws.RIFWorkspace(args['<workspace>'], target)
            workspace.make_dirs()
            pose = pose_from_file(workspace.initial_target_path)
            chain = None
            if chainmap:
                chain = chainmap[workspace.focus_name]
            elif args['--chain']:
                chain = args['--chain']

            if chain:
                print('MAKING PATCHES FOR CHAIN {}'.format(chain))
                poses = []
                for i in range(1, pose.num_chains() + 1):
                    chainpose = pose.split_by_chain(i)
                    info = chainpose.pdb_info().pose2pdb(1)
                    if info.split(' ')[1] in chain and chainpose.residue(1).is_protein():
                        if chainpose.size() < 5:
                            raise('Error: chain {} too small.'.format(chain))
                        else:
                            poses.append(chainpose)
                pose = poses[0]
                if len(poses) > 1:
                    for chainpose in poses[1:]:
                        append_pose_to_pose(pose, chainpose)

            else:
                pose = pose.split_by_chain(1)

            target_pdb = workspace.target_path
            pose.dump_pdb(target_pdb)
            script_path = os.path.join(workspace.patchman_path,
                    'split_to_motifs.py')
            utils.run_command([workspace.python_path,
                workspace.script_path, workspace.target_path])

        except Exception as e:
            print("Error finding patches for {}. Error was:".format(target))
            print(e)
