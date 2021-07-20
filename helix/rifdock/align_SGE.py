import pymol
import sys, os, glob
from helix import workspace as ws


workspace = ws.workspace_from_dir(sys.argv[1])
if 'SGE_TASK_ID' in os.environ:
    task = int(os.environ['SGE_TASK_ID']) - 1
else:
    task = int(sys.argv[2]) - 1
targets = workspace.targets
target = targets[task]
workspace = RIFWorkspace(workspace.root_dir, target)
# folders = sorted(glob.glob(parent + '/*_output'))
# folder = folders[task - 1]
print('Aligning for folder {}'.format(target))
pdbs = sorted(glob.glob(workspace.focus_dir + '/*/docked_full/*.pdb.gz'))

for pdb in pdbs:
    print('Aligning {} to {}'.format(pdb, workspace.target_path))
    pymol.cmd.reinitialize()
    target = pymol.cmd.load(workspace.target_path, 'target')
    mobile = pymol.cmd.load(pdb, 'mobile')
    pymol.cmd.align('mobile and not chain A', 'target')
    pymol.cmd.save(pdb, 'mobile')
