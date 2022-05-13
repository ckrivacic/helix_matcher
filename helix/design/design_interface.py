'''
Design the interface and surroundings of a matched scaffold.
Requires score dataframe and match dataframe.

Usage:
    design_interface.py <workspace> [options]

Options:
    --nstruct=INT, -n  How many designs to make for each input
    --task=INT  Run  a specific task
'''
import sys, os, glob
import pandas as pd
from helix import workspace as ws
from helix import big_jobs
from helix.utils.numeric import euclidean_distance
from helix.utils import utils
from helix.patchman.design_patchman import select_good_residues
import docopt
# Basic PyRosetta stuff
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
# Packing stuff
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.pack.task.residue_selector import ClashBasedShellSelector
# Movers
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.kinematics import MoveMap
# Constraint stuff
from pyrosetta.rosetta.core.scoring.constraints import CoordinateConstraint
from pyrosetta.rosetta.protocols import constraint_generator
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


class InterfaceDesign(object):
    '''Class to help manage things for interface design'''
    def __init__(self, workspace, df_group):
        self.workspace = workspace
        self.df = df_group
        self.design_pose = pose_from_file(
            os.path.join(
                self.workspace.root_dir, self.df.iloc[0].superimposed_file
            )
        )
        self.sfxn_cst = create_score_function('ref2015')
        self.sfxn_cst.set_weight(ScoreType.coordinate_constraint, 1.0)

        self.summarized_residue_scores = utils.safe_load(workspace.find_path(
            'summarized_res_scores.pkl'
        ))

    def setup_task_factory(self):
        '''Set up task factory for design'''
        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())

    def transfer_residue(self, pose2, pose1_resnum, pose2_resnum):
        '''Transfers a rotamer from pose2 to pose1'''
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
        mover = PackRotamersMover(self.sfxn_cst, packertask)
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
        mover = PackRotamersMover(self.sfxn_cst, packertask)
        mover.apply(self.design_pose)


    def prep_design(self):
        '''
        1. Open dataframe, exported scaffold pose, and helix pose
        2. Transfer residues/ROTAMERS over to scaffold pose
            a. For each helix residue, if <0.5 A from scaffold residue, transfer rotamer over with residue type bonus

        Expects a dataframe group where all rows have the same "superimposed_file"
        '''

        final_repack_resis = []
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
                    if closest_row.dist < 0.5:
                        self.transfer_residue(helix_pose, closest_row.res2, helix_index)
                        final_repack_resis.append(int(closest_row.res2))


        # Now that all residues have been transferred, repack and minimize
        # For task factory, repack all the interface residues and a clash-based repack shell. Use constraints.
        self.special_residues = final_repack_resis
        print('Constraints present for the following residues:')
        print(', '.join([str(x) for x in final_repack_resis]))
        idx_selector = utils.list_to_res_selector(final_repack_resis)
        clash_selector = ClashBasedShellSelector()
        clash_selector.set_focus(idx_selector)
        clash_selector.set_include_focus(True)
        clash_selector.set_num_shells(2)
        clash_selector.invert(True)

        tf = TaskFactory()
        no_packing = operation.PreventRepackingRLT()
        static = operation.OperateOnResidueSubset(no_packing,
                                                  clash_selector)
        tf.push_back(static)

        # Repack with constraints
        packertask = tf.create_task_and_apply_taskoperations(self.design_pose)
        mover = PackRotamersMover(self.sfxn_cst, packertask)
        mover.apply(self.design_pose)

        # Add backbone constraints for minimization
        coord_cst = constraint_generator.CoordinateConstraintGenerator()
        # if not self.sc_cst:
        #     coord_cst.set_sidechain(False)
        constraints = coord_cst.apply(self.design_pose)
        for cst in constraints:
            self.design_pose.add_constraint(cst)

        # Minimize with constraints
        movemap = self.setup_movemap()
        minmover = MinMover()
        minmover.movemap(movemap)
        minmover.score_function(self.sfxn_cst)
        minmover.apply(self.design_pose)

        # self.design_pose.dump_pdb('testout.pdb')

    def setup_movemap(self):
        movemap = MoveMap()
        movemap.set_bb(False)
        movemap.set_chi(False)
        for idx, row in self.df.iterrows():
            start = row.rosetta_helix_start
            stop = row.rosetta_helix_stop
            movemap.set_bb_true_range(start, stop)
            movemap.set_chi_true_range(start, stop)

        return movemap


def test_prep():
    df = pd.read_pickle(sys.argv[1])
    init()
    workspace = ws.workspace_from_dir(sys.argv[1])
    for name, group in df.groupby('superimposed_file'):
        print(name)
        setup = InterfaceDesign(workspace, group)
        setup.prep_design()


def main():
    args = docopt.docopt(__doc__)
    try:
        workspace, job_info = big_jobs.initiate()
        job_info['inputs']
    except:
        print('Maybe this is local?', flush=True)
        workspace = ws.workspace_from_dir(args['<workspace>'])
        job_info = {
            'task_id': args['--task'],
        }
    input_df_path = os.path.join(workspace.complex_dir, 'exported_models.pkl')
    input_df = utils.safe_load(input_df_path)
    inputs = sorted(input_df['superimposed_file'].unique())
    if not job_info['task_id']:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0
    else:
        task_id = int(job_info['task_id'])

    dalphaball = os.path.join(workspace.rosetta_dir,
                              'source', 'external', 'DAlpahBall',
                              'DAlphaBall.gcc')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
         ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false'.format(dalphaball))

    nstruct = int(args['--nstruct'])
    total_jobs = len(inputs) * nstruct
    print('TOTAL JOBS: {}'.format(total_jobs), flush=True)

    # Check if we already completed this task, and if so, exit.
    pickle_outdir = workspace.design_dir
    dataframe_out = os.path.join(pickle_outdir, f'task_{task_id}.pkl')

    if os.path.exists(dataframe_out):
        print(f'Task already completed ({dataframe_out} exists). \nExiting.', flush=True)
        sys.exit(0)

    idx = task_id % len(inputs)
    group_name = inputs[idx]
    group = input_df[input_df['superimposed_file'] == group_name]

    designer = InterfaceDesign(workspace, group)

if __name__=='__main__':
    test_prep()