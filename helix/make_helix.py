import sys, os
from ss_generator.alpha_helix import build_ideal_straight_alpha_helix as build
from pyrosetta import *
import numpy as np

'''
Usage:
    make_helix <nres> [outpath]
'''

def make_helix(pose):
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
    phi = -57
    psi = -47
    ca_n_c_o_torsion = 0
    omega = 180

    # Make geometries helical
    for i in range(1, pose.size() + 1):
        pose.set_phi(i, phi)
        pose.set_psi(i, psi)
        pose.set_omega(i, omega)

if __name__=='__main__':

    length = int(sys.argv[1])
    if len(sys.argv) > 2:
        out = os.path.join(
                sys.argv[2], '{}turn_dock_helix.pdb'.format(
                    round(int(length)/3.5))
                )
    else:
        out = './{}turn_dock_helix.pdb'.format(
                round(int(length)/3.5))

    coords = build(length)
    init()
    pose = rosetta.core.pose.Pose()
    typeset = pose.residue_type_set_for_pose()
    new_rsd = rosetta.core.conformation.ResidueFactory.create_residue(
            typeset.name_map("VAL")
            )

    pose.append_residue_by_jump(new_rsd, 1)
    for i in range(1, length):
        pose.conformation().safely_append_polymer_residue_after_seqpos(new_rsd,
                i, True)

    make_helix(pose)
    print('Saving to {}'.format(out))
    pose.dump_pdb(out)
