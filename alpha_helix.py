import numpy as np

from . import geometry
from . import basic


def build_ideal_straight_alpha_helix(length):
    '''Build an ideal straight alpha helix.'''
    
    # Set the basic parameters

    c_n_length = 1.32869
    n_ca_length = 1.458
    ca_c_length = 1.52326
    c_o_length = 1.24
    c_n_ca_angle = np.radians(121.7)
    n_ca_c_angle = np.radians(111.2)
    ca_c_n_angle = np.radians(116.2)
    n_c_o_angle = np.radians(125)
    ca_c_o_angle = np.radians(121)
    phi = np.radians(-57)
    psi = np.radians(-47)
    ca_n_c_o_torsion = np.radians(0)
    omega = np.radians(180)

    # Build the first residue

    helix = [{'ca' : np.array([0, 0, 0]),
              'n' : n_ca_length * np.array([np.sin(n_ca_c_angle), np.cos(n_ca_c_angle), 0]),
              'c' : np.array([0, ca_c_length, 0])}]

    # Build the rest of residues

    for i in range(1, length):
        res = {}
        res['n'] = geometry.cartesian_coord_from_internal_coord(
                    helix[-1]['n'], helix[-1]['ca'], helix[-1]['c'],
                    c_n_length, ca_c_n_angle, psi)
        
        res['ca'] = geometry.cartesian_coord_from_internal_coord(
                    helix[-1]['ca'], helix[-1]['c'], res['n'],
                    n_ca_length, c_n_ca_angle, omega)

        res['c'] = geometry.cartesian_coord_from_internal_coord(
                    helix[-1]['c'], res['n'], res['ca'],
                    ca_c_length, n_ca_c_angle, phi)

        helix.append(res)

    # Add oxygen atoms

    for i in range(len(helix) - 1):
        helix[i]['o'] = geometry.cartesian_coord_from_internal_coord(
                        helix[i + 1]['ca'], helix[i + 1]['n'], helix[i]['c'],
                        c_o_length, n_c_o_angle, ca_n_c_o_torsion)

    # Add oxgen for the last residue

    helix[-1]['o'] = geometry.cartesian_coord_from_internal_coord(
                        helix[-1]['n'], helix[-1]['ca'], helix[-1]['c'],
                        c_o_length, ca_c_o_angle, np.radians(133))

    return helix[:length]

def helix_direction(res1, res2, res3):
    '''Get the helix direction from 3 consecutive residues.'''

    # Get the peptide bond frames

    frame1 = geometry.create_frame_from_three_points(
                res1['c'], res2['n'], res2['ca'])
    frame2 = geometry.create_frame_from_three_points(
                    res2['c'], res3['n'], res3['ca'])

    return geometry.rotation_matrix_to_axis_and_angle(
            np.dot(np.transpose(frame2), frame1))[0]

def get_peptide_bond_rotation_angle(axis, ca_c, n_ca, ref_ca_c):
    '''Get the peptide bond rotation angle given the 
    rotation axis, the ca_c vector and n_ca vector of
    this peptide bond. Choose the solution such that
    the new ca_c direction is near the ref_ca_c direction.
    '''
    axis = geometry.normalize(axis)
    ca_c = geometry.normalize(ca_c)
    n_ca = geometry.normalize(n_ca)
    ref_ca_c = geometry.normalize(ref_ca_c)

    theta1 = np.pi - np.radians(111.2) # The N_CA_C angle
    theta2 = geometry.angle(ca_c, axis)

    # Get the position directions of the new ca_c vector

    ps = geometry.intersections_of_circles_on_unit_sphere(
            n_ca, axis, theta1, theta2)
    new_ca_c = None

    if ps is None:
        print("WARNING! Cannot build along required direction. Use reference direction.")
        new_ca_c = ref_ca_c

    elif np.linalg.norm(ref_ca_c - ps[0]) < np.linalg.norm(ref_ca_c - ps[1]):
        new_ca_c = ps[0]

    else:
        new_ca_c = ps[1]

    # Get the rotation angle

    p1 = geometry.normalize(ca_c - np.dot(ca_c, axis) * axis)
    p2 = geometry.normalize(new_ca_c - np.dot(new_ca_c, axis) * axis)

    return np.sign(np.dot(axis, np.cross(p1, p2))) * np.arccos(np.dot(p1, p2))

def build_alpha_helix_from_directions(directions):
    '''Build an alpha helix from a list of directions.
    The number of residues will be 2 + number of directions.
    '''

    helix = build_ideal_straight_alpha_helix(len(directions) + 2)

    # Align the direction defined by the first 3 residues to the first direction

    d0 = helix_direction(*helix[:3])
    M0 = geometry.rotation_matrix_to_superimpose_two_vectors(d0, directions[0])
    helix = basic.transform_residue_list(helix, M0, np.zeros(3))

    # Change the directions of residues

    for i in range(1, len(directions)):
        res_id = i + 1

        # Get the rotation angle that make the n_ca_c angle ideal

        theta = get_peptide_bond_rotation_angle(directions[i],
                helix[res_id - 1]['c'] - helix[res_id - 1]['ca'],
                helix[res_id]['ca'] - helix[res_id]['n'],
                helix[res_id]['c'] - helix[res_id]['ca'])
    
        # Get the transformation

        frame1 = geometry.create_frame_from_three_points(
                helix[res_id - 1]['ca'], helix[res_id - 1]['c'], helix[res_id]['n'])
        
        frame2 = geometry.create_frame_from_three_points(
                helix[res_id]['ca'], helix[res_id]['c'], helix[res_id + 1]['n'])

        m = geometry.rotation_matrix_from_axis_and_angle(directions[i], theta)
        M = np.dot(np.dot(m, np.transpose(frame1)), frame2)
        
        t = helix[res_id]['ca'] - np.dot(M, helix[res_id]['ca'])

        # Transform the current residue
       
        for atom in ['c', 'o']:
            helix[res_id][atom] = np.dot(M, helix[res_id][atom]) + t

        # Transform the rest of the strand

        for j in range(res_id + 1, len(helix)):
            helix[j] = basic.transform_residue(helix[j], M, t)
        
    return helix
