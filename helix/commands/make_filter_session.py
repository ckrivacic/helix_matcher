'''
Make a PyMOL session where unfiltered helices have a low opacity, filtered ones have higehr opacity.
Meant to show coverage of a protein's surface after filtering.

Usage:
    helix make_filter_session <rif_folder> [options]

'''
import docopt
import pymol
import os
import helix.workspace as ws
import glob
from helix.utils import colors


def color_as_filtered(object_name):
    color = colors.pymol_color(colors.palette['orange'])
    pymol.cmd.color(color, object_name)


def color_as_base(object_name):
    color = colors.pymol_color(colors.palette['red'])
    pymol.cmd.color(color, object_name)
    pymol.cmd.set('cartoon_transparency', 0.7, object_name)


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<rif_folder>'])
    all_files = glob.glob(os.path.join(workspace.focus_dir,
                                       'patch_*', '*', 'docked_full', '*.pdb.gz'))
    filtered_files = glob.glob(os.path.join(workspace.focus_dir, 'filtered', '*.pdb.gz'))
    filtered_filenames = [os.path.basename(x) for x in filtered_files]

    i = 0
    j = 0
    total = len(all_files)
    for f in all_files:
        loadme = os.path.basename(f) in filtered_filenames
        if loadme:
            j += 1
        i += 1
        if loadme and not j%10 == 0:
            continue
        elif not i%50==0 and not loadme:
            continue
        print(f'Loading file {i} of {total}', end='\r')
        basename = os.path.basename(f).split('.')[0]
        pymol.cmd.load(f, basename)
        if os.path.basename(f) in filtered_filenames:
            color_as_filtered(basename)
        else:
            color_as_base(basename)

    pymol.cmd.remove('chain A')

    pymol.cmd.load(workspace.target_path_clean, 'target')
    white = colors.pymol_color(colors.palette['white'])
    pymol.cmd.color(white, 'target')
    pymol.cmd.hide('cartoon', 'target')
    pymol.cmd.show('surface', 'target')

    pymol.cmd.save('filtered_session.pse')