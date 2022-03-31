'''
Design the interface and surroundings of a matched scaffold.
Requires score dataframe and match dataframe.
'''
import sys, os
import pandas as pd
from helix import workspace as ws
from helix.utils.numeric import euclidean_distance
from helix.utils import utils
# Basic PyRosetta stuff
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
# Packing stuff
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
# Constraint stuff
from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring.func import HarmonicFunc
from pyrosetta.rosetta.core.scoring import ScoreType


def get_resmap(pose1, pose2, pose1_start, pose1_stop, pose2_start, pose2_stop):
    '''
    Make a map of CA distances for all residues in pose1 to all residues in pose2
    '''
    resmap = {}
    for res1 in range(pose1_start, pose1_stop + 1):
        if res1 not in resmap:
            resmap[res1] = {}
        for res2 in range(pose2_start, pose2_stop + 1):
            xyz1 = pose1.residue(res1).xyz('CA')
            xyz2 = pose2.residue(res2).xyz('CA')
            resmap[res1][res2] = euclidean_distance(xyz1, xyz2)

    if len(resmap) > 0:
        resmap = pd.DataFrame(resmap).fillna(0).unstack().reset_index()
        resmap.columns = ['res1', 'res2', 'dist']
    else:
        resmap = None

    return resmap


class InterfaceSetup(object):
    '''Class to help manage things for interface design'''
    def __init__(self, workspace, df_group):
        self.workspace = workspace
        self.df = df_group
        self.design_pose = pose_from_file(
            os.path.join(
                self.workspace.root_dir, self.df.iloc[0].superimposed_file
            )
        )
        self.sfxn = create_score_function('ref2015')
        self.sfxn.set_weight(ScoreType.coordinate_constraint, 1.0)

    def transfer_residue(self, pose2, pose1_resnum, pose2_resnum, special_res=False):
        '''Transfers a rotamer from pose2 to pose1'''
        if special_res:
            current_rsd_type_ptr = pose2.residue_type_ptr(position)
            new_rsd_type_mutable = rosetta.core.chemical.MutableResidueType(current_rsd_type_ptr)
            new_rsd_type_mutable.add_variant_type(rosetta.core.chemical.SPECIAL_ROT)
            new_rsd_type = rosetta.core.chemical.ResidueType.make(new_rsd_type_mutable)
            rosetta.core.pose.replace_pose_residue_copying_existing_coordinates(self.design_pose,
                                                                                pose1_resnum, new_rsd_type)
        else:
            pose2_residue = pose2.residue(pose2_resnum)

            # Define operations
            idx_selector = residue_selector.ResidueIndexSelector(int(pose1_resnum))
            restypes = operation.RestrictToSpecifiedBaseResidueTypes(
                utils.strlist_to_vector1_str([pose2_residue.name3()]),
                idx_selector)
            not_selector = residue_selector.NotResidueSelector(idx_selector)
            no_packing = operation.PreventRepackingRLT()
            static = operation.OperateOnResidueSubset(no_packing,
                                                      not_selector)
            # Apply operations to task factory
            tf = TaskFactory()
            tf.push_back(operation.InitializeFromCommandline())
            tf.push_back(restypes)
            tf.push_back(static)
            packertask = tf.create_task_and_apply_taskoperations(self.design_pose)

            # Pack to mutate residue
            mover = PackRotamersMover(self.sfxn, packertask)
            mover.apply(self.design_pose)

            # Now add constraints
            for atom in range(pose2_residue.first_sidechain_atom(), pose2_residue.natoms()):
                if pose2_residue.atom_is_hydrogen(atom):
                    continue
                xyz = pose2_residue.xyz(atom)
                id = AtomID(atom, int(pose1_resnum))
                reference = AtomID(1, 1)
                func = HarmonicFunc(0, 1)
                cst = CoordinateConstraint(id, reference, xyz, func)
                self.design_pose.add_constraint(cst)

            # Repack with constraints
            packertask = tf.create_task_and_apply_taskoperations(self.design_pose)
            mover = PackRotamersMover(self.sfxn, packertask)
            mover.apply(self.design_pose)


    def prep_design(self):
        '''
        1. Open dataframe, exported scaffold pose, and helix pose
        2. Transfer residues/ROTAMERS over to scaffold pose
            a. For each helix residue, if <0.5 A from scaffold residue, transfer rotamer over with residue type bonus

        Expects a dataframe group where all rows have the same "superimposed_file"
        '''

        for idx, row in self.df.iterrows():
            # Helix residue positions are for just that chain
            docked_pose = pose_from_file(
                os.path.join(
                    self.workspace.root_dir, row.design_file
                )
            )
            helix_pose = docked_pose.split_by_chain(2)

            resmap = get_resmap(helix_pose, self.design_pose, row.start, row.stop, row.rosetta_helix_start, row.rosetta_helix_stop)

            # print(row.interfac_residues)
            # Interface residues, on the other hand, are numbered by rosetta number of the whole docked pose, so we need
            # to adjust those #s
            for resi in row.interfac_residues:
                if docked_pose.chain(resi) == 2:
                    helix_index = resi - docked_pose.chain_begin(2) + 1
                    if (helix_index > row.stop or helix_index < row.start):
                        continue
                    res_distances = resmap[resmap.res1 == helix_index]
                    # print(res_distances)
                    closest_row = res_distances.sort_values(by='dist').iloc[0]
                    # print(closest_row.res2)
                    # print(closest_row.dist)
                    if closest_row.dist < 1.0:
                        self.transfer_residue(helix_pose, closest_row.res2, helix_index)

        self.design_pose.dump_pdb('testout.pdb')



def test_prep():
    df = pd.read_pickle(sys.argv[1])
    init()
    workspace = ws.workspace_from_dir(sys.argv[1])
    for name, group in df.groupby('superimposed_file'):
        print(name)
        setup = InterfaceSetup(workspace, group)
        setup.prep_design()

if __name__=='__main__':
    test_prep()