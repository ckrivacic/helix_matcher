import os, glob, sys
from shutil import copyfile

lengths = [3,4,6,8]
folders = glob.glob(sys.argv[1]) + '/*/cluster_representatives/'

for f in folders:
    for length in lengths:
        subfolder = '{}_turn'.format(length)
        for pdb in glob.glob(
                os.path.join(
                    f, subfolder
                    )
                ):
            os.symlink(pdb, 
                    os.path.join(
                        f, os.path.basename(pdb)
                        ))
