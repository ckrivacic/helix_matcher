import pymol
import sys, os, glob


parent = sys.argv[1]
for folder in glob.glob(parent + '/*_output'):
    print('Aligning for folder {}'.format(folder))
    pdbs = sorted(glob.glob(folder + '/*/docked_full/*.pdb.gz'))

    for pdb in pdbs:
        print('Aligning {} to {}'.format(pdb, pdbs[0]))
        pymol.cmd.reinitialize()
        target = pymol.cmd.load(pdbs[0], 'target')
        mobile = pymol.cmd.load(pdb, 'mobile')
        pymol.cmd.align('mobile and chain B', 'target and chain B')
        pymol.cmd.save(pdb, 'mobile')
