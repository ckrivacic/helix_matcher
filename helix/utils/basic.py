import numpy as np

from helix.utils import geometry


def transform_residue(residue, M, t):
    '''Transform a residue by a rotation M and a
    translation t. Return the transformed residue.
    '''
    new_res = {}

    for key in residue.keys():
        new_res[key] = np.dot(M, residue[key]) + t

    return new_res

def transform_residue_list(res_list, M, t):
    '''Transform a residue list by a rotation M and a
    translation t. Return the new list.
    '''
    return [transform_residue(res, M, t) for res in res_list]

def get_phi(strand, res_id):
    '''Get the phi torsions of a residue.'''
    return geometry.dihedral(strand[res_id - 1]['c'], strand[res_id]['n'],
            strand[res_id]['ca'], strand[res_id]['c'])

def get_psi(strand, res_id):
    '''Get the psi torsions of a residue.'''
    return geometry.dihedral(strand[res_id]['n'], strand[res_id]['ca'],
            strand[res_id]['c'], strand[res_id + 1]['n'])

def change_torsions(strand, res_id, phi, psi):
    '''Change the phi, psi angles of a residue in
    a strand. The input torsions should be in radians.
    '''

    # Rotate the psi torsion

    if 0 <= res_id < len(strand) - 1:
    
        psi_old = geometry.dihedral(strand[res_id]['n'], strand[res_id]['ca'], 
                strand[res_id]['c'], strand[res_id + 1]['n'])

        # Get the rotation matrix

        axis = strand[res_id]['c'] - strand[res_id]['ca']
        M = geometry.rotation_matrix_from_axis_and_angle(axis, psi - psi_old)
        t = strand[res_id]['ca'] - np.dot(M, strand[res_id]['ca'])

        # Rotate subsequent atoms

        strand[res_id]['o'] = np.dot(M, strand[res_id]['o']) + t

        for i in range(res_id + 1, len(strand)):
            strand[i] = transform_residue(strand[i], M, t)

    # Rotate the phi torsion
    
    if 0 < res_id < len(strand):

        phi_old = geometry.dihedral(strand[res_id - 1]['c'], strand[res_id]['n'],
                strand[res_id]['ca'], strand[res_id]['c'])

        # Get the rotation matrix

        axis = strand[res_id]['ca'] - strand[res_id]['n']
        M = geometry.rotation_matrix_from_axis_and_angle(axis, phi - phi_old)
        t = strand[res_id]['ca'] - np.dot(M, strand[res_id]['ca'])

        # Rotate subsequent residues

        for key in strand[res_id].keys():
            if key != 'h':
                strand[res_id][key] = np.dot(M, strand[res_id][key]) + t

        for i in range(res_id + 1, len(strand)):
            strand[i] = transform_residue(strand[i], M, t)

def get_hb_co_coord_from_nh(n, h):
    '''Get the ideal coordinates of C and O from
    coordinates of N and H.
    '''
    v = geometry.normalize(h - n)
    return (n + 3.9 * v, n + 2.7 * v)

def get_hb_nh_coord_from_co(c, o):
    '''Get the ideal coordinates of N and H from
    coordinates of C and O.
    '''
    v = geometry.normalize(o - c)
    return (c + 3.9 * v, c + 2.9 * v)

def get_peptide_bond_transformation(phi, psi):
    '''Get the rotation matrix and translation vector of a peptide bond
    corresponding to a pair of phi psi torsions. The reference frame
    is build on the C, N and CA atoms.
    '''
    # Set the parameters

    n_ca_length = 1.47
    ca_c_length = 1.53
    c_n_length = 1.32
    n_ca_c_angle = np.radians(111.2)
    ca_c_n_angle = np.radians(114)
    c_n_ca_angle = np.radians(123)
    omega = np.pi

    # Get the coordinates

    c1 = c_n_length * np.array([np.sin(c_n_ca_angle), np.cos(c_n_ca_angle), 0])
    n1 = np.array([0, 0, 0])
    ca1 = np.array([0, n_ca_length, 0])
    c2 = geometry.cartesian_coord_from_internal_coord(c1, n1, ca1, 
            ca_c_length, n_ca_c_angle, phi)
    n2 = geometry.cartesian_coord_from_internal_coord(n1, ca1, c2,
            c_n_length, ca_c_n_angle, psi)
    ca2 = geometry.cartesian_coord_from_internal_coord(ca1, c2, n2,
            n_ca_length, c_n_ca_angle, omega)

    # Get the transformation

    return np.transpose(geometry.create_frame_from_three_points(c2, n2, ca2)), n2 - n1
