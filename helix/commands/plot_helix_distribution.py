#!/usr/bin/env python3
'''Plot the distribution of helices in the beta sheet coordinates.
Adapted from LUCS, by Xingjie Pan.

Usage:
    helix plot_helix_distribution <workspace> [options]

Options:
    --target=STR  Only plot for the given target
    --overwrite, -o  Recalculate
    --task=INT  Task number
    --noplot  Don't actually plot
    --load  Load
    --save  Save the plot
'''

import os
import json
import glob

import docopt
import helix.workspace as ws
from helix.utils.colors import palette

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches

import pyrosetta
from pyrosetta import rosetta


def xyz_to_np_array(xyz):
    '''Convert an xyz vector to a numpy array.'''
    return np.array([xyz.x, xyz.y, xyz.z])


def get_sheet_ca_positions(pose, sheet_residues):
    '''Get CA positions of a beta sheet'''
    sheet_ca_positions = []

    for strand in sheet_residues:
        strand_ca_positions = []

        for seqpos in strand:
            if not (seqpos is None):
                xyz = pose.residue(seqpos).xyz('CA')
                strand_ca_positions.append(xyz_to_np_array(xyz))
            else:
                strand_ca_positions.append(None)

        sheet_ca_positions.append(strand_ca_positions)

    return sheet_ca_positions


def get_sheet_residue_local_frames(pose, sheet_residues, strand_directions, z_ref_residue):
    '''Get local frames for each beta sheet residues.
    The x direction is determined by the strand direction.
    The y direction is determined by neighbor strands.
    The z direction is determined by positions of CA and CB atoms.
    '''
    sheet_res_frames = []

    for i, strand in enumerate(sheet_residues):
        strand_res_frames = []

        for j, seqpos in enumerate(strand):
            if not (seqpos is None):
                # Get the Z direction

                ca_xyz = pose.residue(seqpos).xyz('CA')

                if pose.residue(seqpos).name3() == 'GLY':
                    cb_xyz = pose.residue(seqpos).xyz('2HA')
                else:
                    cb_xyz = pose.residue(seqpos).xyz('CB')

                z_axis_no_sign = xyz_to_np_array((cb_xyz - ca_xyz).normalized())
                z_sign = 1 if (0 == (j - z_ref_residue) % 2) else -1
                z_axis = z_sign * z_axis_no_sign

                # Get the Y direction

                y_axis_no_sign = xyz_to_np_array((pose.residue(seqpos).xyz('C')
                                                  - pose.residue(seqpos).xyz('N')).normalized())

                y_sign = 1 if strand_directions[i] else -1
                y_axis_before_correction = y_sign * y_axis_no_sign
                y_axis_before_correction = y_axis_before_correction - np.dot(y_axis_before_correction, z_axis) * z_axis

                y_axis = y_axis_before_correction / np.linalg.norm(y_axis_before_correction)

                # Get the X direction

                x_axis = np.cross(y_axis, z_axis)

                strand_res_frames.append([x_axis, y_axis, z_axis])

            else:
                strand_res_frames.append(None)

        sheet_res_frames.append(strand_res_frames)

    return sheet_res_frames


def project_point_to_sheet(sheet_ca_positions, sheet_res_frames, point, vector=None):
    '''Project a point onto the sheet coordinates.
    Use a gaussian function to weight the contributations
    of different local frames.
    Also project a vector attached to that point if given.
    '''
    ref_vertical_length = 3.3
    ref_horizontal_length = 4.6

    ca_p_linear = [ca for strand in sheet_ca_positions for ca in strand if not (ca is None)]

    res_frames_linear = [np.array(f) for strand in sheet_res_frames for f in strand if not (f is None)]

    sheet_coords_linear = [np.array([i * ref_horizontal_length, j * ref_vertical_length])
                           for i in range(len(sheet_ca_positions)) for j in range(len(sheet_ca_positions[0]))
                           if not (sheet_ca_positions[i][j] is None)]

    # Assign weights to sheet residues

    dists = [np.linalg.norm(point - ca) for ca in ca_p_linear]
    dist_scale = 5 ** 2
    weights = [np.exp(-d * d / dist_scale) for d in dists]
    total_weight = sum(weights)
    weights = [w / total_weight for w in weights]

    # Get the average positions

    mean_position = sum(weights[i] * ca_p_linear[i] for i in range(len(ca_p_linear)))
    mean_sheet_coord = sum(weights[i] * sheet_coords_linear[i] for i in range(len(sheet_coords_linear)))
    mean_frame = sum(weights[i] * res_frames_linear[i] for i in range(len(res_frames_linear)))

    # Get the projected sheet coordinates

    diff = point - mean_position

    x = mean_sheet_coord[0] + np.dot(diff, mean_frame[0])
    y = mean_sheet_coord[1] + np.dot(diff, mean_frame[1])

    # Project the attached vector

    if not (vector is None):
        v_x = np.dot(vector, mean_frame[0])
        v_y = np.dot(vector, mean_frame[1])
        return np.array([x, y]), np.array([v_x, v_y])

    return np.array([x, y])


def project_a_helix_to_sheet_coords(pose, helix, sheet_ca_positions, sheet_res_frames):
    '''Project A helix to sheet coordinates.
    The helix is defined as a pair (start, stop).
    Return the center of helix and the direction of the helix
    in the sheet coordinates.
    '''
    helix_cas = [xyz_to_np_array(pose.residue(i).xyz('CA')) for i in range(helix[0], helix[1] + 1)]
    helix_dirs = [xyz_to_np_array(pose.residue(i).xyz('O') - pose.residue(i).xyz('C'))
                  for i in range(helix[0], min(helix[1] + 1, pose.size()))]

    center = sum(helix_cas) / len(helix_cas)
    center_dir = sum(helix_dirs)

    c_pj, d_pj = project_point_to_sheet(sheet_ca_positions, sheet_res_frames, center, center_dir)

    return c_pj, d_pj / np.linalg.norm(d_pj)


def plot_helices(pose, helix_coords, sheet_ca_positions, sheet_res_frames, match=False):
    '''Plot helices on under the sheet coordinates.
    A helix coordinate is defined by 2 2D vectors
    (helix_center_position, helix_direction)
    '''
    X = []
    Y = []
    U = []
    V = []

    for hc in helix_coords:
        c_pj, d_pj = hc
        X.append(c_pj[0] - d_pj[0] / 2)
        Y.append(c_pj[1] - d_pj[1] / 2)
        U.append(d_pj[0])
        V.append(d_pj[1])

    color = palette['blue']
    width = 0.002
    headaxislength = 4.5
    headwidth = 10
    headlength = 6
    scale = 30
    alpha = 0.2
    lw = 0.2
    ec = palette['white']
    if match:
        # headwidth = 3
        # headlength = 3
        # headaxislength = 2
        # scale = 30
        width = 0.002
        color = palette['teal']
        alpha = 1
        # lw=0.2
        ec = palette['white']
    # else:
    #     color = palette['blue']
    #     width = 0.001
    #     headaxislength=4.5
    #     headwidth = 20
    #     headlength=12
    #     scale = 30
    #     alpha = 0.2
    #     lw=0.1
    #     ec=palette['blue']
    f = plt.figure()
    ax = plt.subplot(111)
    plt.quiver(X, Y, U, V, width=width, color=color, headwidth=headwidth, headlength=headlength,
               headaxislength=headaxislength, scale=scale, alpha=alpha, lw=lw, ec=ec)
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")

    # plt.show()


def plot_the_underlying_sheet(sheet_ca_positions, sheet_res_frames):
    '''Plot the underlying sheet.'''
    sheet_coords = []

    # Get the projected coordinates of sheet residues

    polygon_color = np.array([1, 1, 1]) * 0.9

    for i in range(len(sheet_ca_positions)):
        strand_coords = []
        for j in range(len(sheet_ca_positions[i])):

            if sheet_ca_positions[i][j] is None:
                continue

            color = polygon_color if j % 2 else 'white'

            strand_coords.append(
                (project_point_to_sheet(sheet_ca_positions, sheet_res_frames, sheet_ca_positions[i][j]),
                 color))

            sheet_coords.append(strand_coords)

    # Plot the poligons of peptide bonds

    ax = plt.axes()

    for i in range(len(sheet_coords)):
        for j in range(len(sheet_coords[i]) - 1):
            l = np.array([-1, 0])
            r = np.array([1, 0])

            p1 = sheet_coords[i][j][0]
            p2 = sheet_coords[i][j + 1][0]

            polygon = plt.Polygon([p1 + l, p1 + r, p2 + r, p2 + l],
                                  edgecolor=polygon_color, facecolor=sheet_coords[i][j][1], zorder=-100)
            ax.add_patch(polygon)

    # plt.show()


def plot_test(sheet_ca_positions, sheet_res_frames, points):
    X = []
    Y = []

    for p in points:
        sheet_coord = project_point_to_sheet(sheet_ca_positions, sheet_res_frames, p)
        X.append(sheet_coord[0])
        Y.append(sheet_coord[1])

    plt.plot(X, Y)


def get_all_helix_coords_for_one_design(pdb_file, insertion_file, sheet_ca_positions, sheet_res_frames):
    '''Get all coordinates for remodeled helices for one design.'''
    pose = rosetta.core.import_pose.pose_from_file(pdb_file)
    dssp_str = rosetta.core.scoring.dssp.Dssp(pose).get_dssp_secstruct()

    # Find all helices
    with open(insertion_file, 'r') as f:
        insertion_info = json.load(f)

    helix_coords = []

    for d in insertion_info:
        h_start = d['stop']
        h_stop = d['start']

        # if h_start > h_stop: continue

        for i in range(d['start'], d['stop'] + 1):
            if dssp_str[i - 1] == 'H':
                if i < h_start:
                    h_start = i
                if i > h_stop:
                    h_stop = i
        if h_start > h_stop: continue

        helix_coords.append(
            project_a_helix_to_sheet_coords(pose, (h_start, h_stop), sheet_ca_positions, sheet_res_frames))

    return helix_coords


def get_json(workspace, pdb_path):
    model_no = os.path.basename(pdb_path).split('_')[1]
    lhl_folder = os.path.join(workspace.root_dir, '..', 'regenerated_data_sets_2020_03',
                              'sequence_design_for_LHL_reshaping_2lv8_two_LHL',
                              'selected_designs_for_state_count')
    insertion_file = os.path.join(lhl_folder, f'insertion_points_{model_no}.json')
    with open(insertion_file, 'r') as f:
        insertions = json.load(f)
    return insertions


def get_all_helix_coords_for_data_set(data_path, sheet_ca_positions, sheet_res_frames, match_filenames, task=None):
    '''Get all coordinates for remodeled helices in a data set.'''
    # Get all pdb files and insertion files
    print(f"TASK {task}")

    pdb_files = []
    insertion_files = []

    for f in os.listdir(data_path):
        if f.endswith('.pdb.gz'):
            pdb_id = f.split('.')[0].split('_')[-1]
            insertion_file = os.path.join(data_path, 'insertion_points_{0}.json'.format(pdb_id))

            if os.path.exists(insertion_file):
                insertion_files.append(insertion_file)
                pdb_files.append(os.path.join(data_path, f))

    pdb_files = sorted(pdb_files)
    # Get the helix coords

    helix_coords = {}
    total = len(pdb_files)
    if task is None:
        task = 0
        increment = 1
    else:
        increment = 500

    for i in range(task, total, increment):
        print(i)
        try:
            these_coords = get_all_helix_coords_for_one_design(pdb_files[i], insertion_files[i], sheet_ca_positions,
                                                                sheet_res_frames)
            helix_coords[os.path.basename(pdb_files[i])] = these_coords
        except:
            print(f'Could not get coords for {pdb_files[i]}')
            continue

    return helix_coords


def dump_helix_coords(helix_coords, file_name):
    '''Dump helix coordinates to a json file'''
    h_coords_serial = {}# [[list(c), list(d)] for c, d in helix_coords]
    for filename in helix_coords:
        h_coords_serial[filename ] = [[list(c), list(d)] for c, d in helix_coords[filename]]

    with open(file_name, 'w') as f:
        json.dump(h_coords_serial, f)


def load_helix_coords(file_name):
    '''Load helix coordinates from a json file'''
    with open(file_name, 'r') as f:
        h_coords_serial = json.load(f)

    helix_coords = {}
    for filename in h_coords_serial:
        helix_coords[filename] = [[np.array(c), np.array(d)] for c, d in h_coords_serial[filename]]

    return helix_coords


# if __name__ == '__main__':
def main():
    pyrosetta.init()

    args = docopt.docopt(__doc__)
    if not args['--noplot']:
        mpl.use('tkagg')
    if args['--task']:
        taskno = int(args['--task'])
    else:
        taskno = None
    workspace = ws.workspace_from_dir(args['<workspace>'])
    ref_pdb = os.path.expanduser('~/intelligent_design/2lv8.pdb')
    if args['--target']:
        targets = [workspace.target_match_path(args['--target'])]
    else:
        targets = workspace.all_match_workspaces

    pose = rosetta.core.import_pose.pose_from_file(ref_pdb)

    if 'SGE_TASK_ID' in os.environ:
        taskno = int(os.environ['SGE_TASK_ID']) - 1

    # Define residues in a beta sheet. The residues should be aligned.
    # Residues that cannot be aligned should be set to None.

    sheet_residues = [
        [None, None, 27, 28, 29, 30, 31],
        [1, 2, 3, 4, 5, 6, 7],
        [51, 52, 53, 54, 55, 56, 57],
        [None, None, 77, 78, 79, 80, 81],
    ]
    strand_directions = [True, True, True, True]

    # The residue on a strand that determines the Z-direction
    # This value must be set correctly such that as the indices
    # increase, the projected x, y coordinates also increase.
    z_ref_residue = 27

    sheet_ca_positions = get_sheet_ca_positions(pose, sheet_residues)
    sheet_res_frames = get_sheet_residue_local_frames(pose, sheet_residues, strand_directions, z_ref_residue)

    if not args['--noplot']:
        plot_the_underlying_sheet(sheet_ca_positions, sheet_res_frames)

    # ca_points = [xyz_to_np_array(pose.residue(i).xyz('CA')) for i in range(1, pose.size() + 1)]
    # plot_test(sheet_ca_positions, sheet_res_frames, ca_points)

    # helices = [(36, 52), (71, 81)]
    # helix_coords = [project_a_helix_to_sheet_coords(pose, helix, sheet_ca_positions, sheet_res_frames)
    #        for helix in helices]

    lhl_folder = os.path.join(workspace.root_dir, '..', 'regenerated_data_sets_2020_03',
                              'sequence_design_for_LHL_reshaping_2lv8_two_LHL',
                              'selected_designs_for_state_count')

    model_names = []
    for target in targets:
        match_workspace = ws.MatchWorkspace(workspace.root_dir, target)
        data_path = match_workspace.complex_dir
        model_names += ['_'.join(os.path.basename(x).split('.')[0].split('_')[:-1]) for x in glob.glob(os.path.join(data_path, '*.pdb.gz'))]
    #     helix_coords = get_all_helix_coords_for_data_set(data_path, sheet_ca_positions, sheet_res_frames)

    print('Checking for path')
    base_helix_coords_path = os.path.join(workspace.root_dir, 'lucs_helix_info', f'all_helices_{taskno}.json')
    if not os.path.exists(os.path.join(workspace.root_dir, 'lucs_helix_info')):
        print('Path does not exist; creating')
        os.makedirs(os.path.join(workspace.root_dir, 'lucs_helix_info'), exist_ok=True)
    # match_helix_coords_path = os.path.join(workspace.root_dir, 'matched_helices.json')
    if os.path.exists(base_helix_coords_path) and not args['--overwrite']:
        print('Loading helices')
        helix_coords = load_helix_coords(base_helix_coords_path)
        # match_coords = load_helix_coords(match_helix_coords_path)
    elif args['--load']:
        print("Loading all helices")
        coord_paths = glob.glob(os.path.join(workspace.root_dir, 'lucs_helix_info', 'all_helices*json'))
        helix_coords = {}
        for path in coord_paths:
            these_coords = load_helix_coords(path)
            helix_coords = {**helix_coords, **these_coords}
    else:
        print(f'Calculating coords; task {taskno}')
        helix_coords = get_all_helix_coords_for_data_set(lhl_folder, sheet_ca_positions, sheet_res_frames,
                                                                       model_names, task=taskno)
        dump_helix_coords(helix_coords, base_helix_coords_path)
        # dump_helix_coords(match_coords, match_helix_coords_path)

    # dump_helix_coords(helix_coords, 'test.json')
    # helix_coords = load_helix_coords('test.json')
    unmatched = []
    matched = []
    for filename in helix_coords:
        if filename.split('.')[0] in model_names:
            matched.extend(helix_coords[filename])
        else:
            unmatched.extend(helix_coords[filename])

    if not args['--noplot']:
        unmatched.extend(matched)
        plot_helices(pose, unmatched, sheet_ca_positions, sheet_res_frames, match=False)
        # plot_helices(pose, matched, sheet_ca_positions, sheet_res_frames, match=True)

        plt.axes().set_aspect('equal')
        if args['--save']:
            if not args['--target']:
                target = 'all'
            plt.savefig(f"/Volumes/GoogleDrive-109095367122122187045/My Drive/Kortemme_Lab/helix/manuscript/supplemental/{os.path.basename(target)}_unmat.png",
                        transparent=True, dpi=600)
        else:
            plt.show()

if __name__=='__main__':
    main()