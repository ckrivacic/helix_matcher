'''
Analyze how many faces each helix ended up binding to.

Usage:
    design_face_analysis.py <workspace> [options]

Options:
'''
import docopt
import helix.workspace as ws
from helix.commands.plot_helix_distribution import get_json
import prody
import numpy as np
import pandas as pd
from helix.utils.numeric import dihedral


def match_faces(interface_resis, insertion):
    '''Helical faces have a pattern of i, i+3, i+4, i+7.
    Get list of helices that are in the same interface as each other?'''
    start = insertion['start']
    stop = insertion['stop']

    started = False
    for resi in sorted(interface_resis):
        if resi > start and resi < stop:
            if not started:
                i = resi
                started = True  # we are now in the helix/insertion
            else:
                diff = resi - i


def is_same_face(ca_cb_1, ca_cb_2):
    '''Check if two resis are on the same face'''


def get_same_face(atoms, insertion, interface_residues):
    selection = atoms.select(f"(name CA or name CB) and resnum {insertion['start']}:{insertion['stop']}")
    hv = selection.getHierView()
    dihedral_map = {}
    start = insertion['start']
    stop = insertion['stop']
    for i, residue in enumerate(hv.iterResidues()):
        dihedral_map[residue.getResnum()] = {}
    for i, residue in enumerate(hv.iterResidues()):
        coords = residue.getCoords()
        if len(coords) < 2:
            continue
        for j, residue_j in enumerate(hv.iterResidues()):
            if i in dihedral_map[residue_j.getResnum()]:
                dihedral_map[residue.getResnum()][residue_j.getResnum()] = dihedral_map[residue_j.getResnum()][residue.getResnum()]
            if i == j:
                dihedral_map[residue.getResnum()][residue_j.getResnum()] = np.nan
                continue
            coords_j = residue_j.getCoords()
            if len(coords_j) < 2:
                dihedral_map[residue.getResnum()][residue_j.getResnum()] = np.nan
                continue
            dihedral_map[residue.getResnum()][residue_j.getResnum()] = dihedral(coords[0], coords[1], coords_j[0], coords_j[1])

    dihedral_map = pd.DataFrame.from_dict(dihedral_map)
    different_faces = 0
    for i, residue in enumerate(hv.iterResidues()):
        resi = residue.getResnum()
        if resi in interface_residues and resi < stop:
            dihedral_angle = dihedral_map.loc[start, resi]
            if dihedral_angle > 90 or dihedral_angle < -90:
                different_faces += 1
    return different_faces


def main():
    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<workspace>'])
    targets = workspace.targets

    for target in targets:
        match_workspace = ws.MatchWorkspace(workspace.root_dir, target)
        scores = match_workspace.get_scores()
        for idx, row in scores.iterrows():
            interface = row['interfac_residues']
            insertions = get_json(workspace, row['design_file'])
            atoms = prody.parsePDB(row['design_file'])
            num_different = get_same_face(atoms, insertions[0], interface)
            print(num_different)


if __name__=='__main__':
    main()