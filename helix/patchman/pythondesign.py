from pyrosetta import rosetta
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.protocols import constraint_generator
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta import *
from helix.utils import utils
import os


class SpecialRotDesign(object):
    def __init__(self, sfxn=None, script=None, nrepeats=5, ramp_cst=False,
                 special_rotamers=[], special_rot_weight=-1.5, taskfactory=None,
                 movemap=None, ramp_down_constraints=False, bb_cst=False, sc_cst=False):
        '''Initialize'''
        self.nrepeats = nrepeats
        if not script:
            self.script = self.get_default_script()
        else:
            self.script = script
        self.script_size = len(self.script)
        if not sfxn:
            self.sfxn = create_score_function('ref2015')
            self.sfxn.set_weight(ScoreType.coordinate_constraint, 1.0)
        else:
            self.sfxn = sfxn

        self.special_rot_weight = special_rot_weight
        self.special_rotamers = special_rotamers
        if len(self.special_rotamers) > 0:
            self.sfxn.set_weight(rosetta.core.scoring.special_rot, self.special_rot_weight)

        if taskfactory:
            self.taskfactory = taskfactory
        else:
            # self.taskfactory = self.setup_default_taskfactory()
            self.taskfactory = None

        if movemap:
            self.movemap = movemap
        else:
            self.movemap = None

        self.ramp_down_constraints = ramp_down_constraints

        self.bb_cst = bb_cst
        self.sc_cst = sc_cst

        self.best_score = 9999

    def set_scorefxn(self, sfxn):
        self.sfxn = sfxn

    def set_movemap(self, movemap):
        self.movemap = movemap

    def ramp_down_constraints(self, ramp):
        '''boolean'''
        self.ramp_down_constraints = ramp

    def set_task_factory(self, tf):
        self.taskfactory = tf

    def setup_default_taskfactory(self):
        '''Set up a default task factory for interface design'''
        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        selector = residue_selector.ChainSelector('B')
        interface_selector_str = \
            '''
        <InterfaceByVector name="interface_selector" nearby_atom_cut="6.5">
           <Chain chains='A'/>
           <Chain chains='B'/>
        </InterfaceByVector>
        '''
        interface_selector = XmlObjects.static_get_residue_selector(interface_selector_str)
        or_selector = residue_selector.OrResidueSelector(selector,
                                                         interface_selector)
        clash_selector = rosetta.core.pack.task.residue_selector.ClashBasedShellSelector(or_selector)
        clash_selector.invert(False)
        clash_selector.set_include_focus(True)
        not_selector = residue_selector.NotResidueSelector(clash_selector)
        no_packing = operation.PreventRepackingRLT()
        no_design = operation.RestrictToRepackingRLT()
        upweight = \
            '''
        <ProteinLigandInterfaceUpweighter name="upweight_interface" interface_weight="1.5" />
        '''
        upweight_taskop = XmlObjects.static_get_task_operation(upweight)
        static = operation.OperateOnResidueSubset(no_packing,
                                                  not_selector)
        notaa = operation.ProhibitSpecifiedBaseResidueTypes(
            utils.strlist_to_vector1_str(['GLY']),
            selector)
        and_selector = residue_selector.AndResidueSelector(interface_selector,
                                                           selector)
        not_interface = residue_selector.NotResidueSelector(and_selector)
        notdesign = operation.OperateOnResidueSubset(no_design,
                                                     not_interface)
        print('Pushing taskops back', flush=True)
        tf.push_back(notaa)
        tf.push_back(notdesign)
        tf.push_back(static)
        tf.push_back(upweight_taskop)
        print('Finished setting up task factory', flush=True)

        return tf

    def setup_default_movemap(self, pose):
        movemap = MoveMap()
        movemap.set_bb_true_range(pose.chain_begin(2),
                                  pose.chain_end(2))
        movemap.set_chi_true_range(pose.chain_begin(2),
                                   pose.chain_end(2))
        return movemap

    def get_default_script(self):
        '''Try to get Rosetta dir from environment'''
        if 'ROSETTA' in os.environ:
            rosettadir = os.environ['ROSETTA']
        else:
            rosettadir = input('Please provide Rosetta main directory '\
                    '(in the future you can set the environmental '\
                    'variable $ROSETTA)\n')
        default_script_path = os.path.join(rosettadir, 'database',
                'sampling', 'relax_scripts', 'MonomerDesign2019.txt')

        return default_script_path

    def parse_script(self, script_path):
        '''Parse a Relax script; return a list of movers?'''
        script = open(script_path, 'r')
        script_steps = []
        for line in script:
            line = line.strip()
            if line.startswith('repeat'):
                continue
            elif line.startswith('coord_cst_weight'):
                script_steps.append(('coord_cst_weight', float(line.split(' ')[-1])))
            elif  line.startswith('scale'):
                script_steps.append((line.split(' ')[0],
                    float(line.split(' ')[-1])
                    ))
            elif line.startswith('min'):
                script_steps.append(('min', float(line.split(' ')[1])))
            elif line.startswith('repack'):
                script_steps.append(('repack',))
            elif line.startswith('accept_to_best'):
                script_steps.append(('accept_to_best',))

        self.script_steps = script_steps

    def pack_rotamers(self, pose):
        packertask = self.taskfactory.create_task_and_apply_taskoperations(pose)
        local_sfxn = self.sfxn.clone()
        print('Packertask created')
        print(packertask)
        local_sfxn(pose)
        local_sfxn.setup_for_packing(pose, packertask.repacking_residues(), packertask.designing_residues())
        packer_neighbor_graph = rosetta.core.pack.create_packer_graph(
            pose, local_sfxn, packertask
        )

        rotamer_sets = rosetta.core.pack.rotamer_set.RotamerSetsFactory.create_rotamer_sets(pose)
        rotamer_sets.set_task(packertask)
        rotamer_sets.initialize_pose_for_rotsets_creation(pose)
        rotamer_sets.build_rotamers(pose, local_sfxn, packer_neighbor_graph)

        for position in self.specialrot_dict:
            if packertask.design_residue(position):
                position_rotamer_set = rotamer_sets.rotamer_set_for_residue(position)

                # Add special rotamer to the appropriate rotamer set
                if int(position_rotamer_set.resid()) == position:
                    print(f"Adding the following restype to position {position}:", flush=True)
                    print(self.specialrot_dict[position].annotated_name(True), flush=True)
                    print("Originally had {} rotamers".format(
                        position_rotamer_set.num_rotamers()
                    ))
                    position_rotamer_set.add_rotamer_into_existing_group(self.specialrot_dict[position])
                    print("Position now has {} rotamers".format(
                        position_rotamer_set.num_rotamers()
                    ))

        local_sfxn.setup_for_packing_with_rotsets(pose, rotamer_sets)
        rotamer_sets.prepare_sets_for_packing(pose, local_sfxn)
        ig = rosetta.core.pack.interaction_graph.InteractionGraphFactory.create_and_initialize_annealing_graph(
            packertask, rotamer_sets, pose, local_sfxn, packer_neighbor_graph)
        rosetta.core.pack.pack_rotamers_run(pose, packertask, rotamer_sets, ig)
        ig.clean_up_after_packing(pose)
        local_sfxn(pose)

        test = False
        if test:
            print(f'Scorefxn for pose: {local_sfxn(pose)}')
            ref = create_score_function('ref2015')
            print(f"Non-special sfxn for pose: {ref(pose)}")
            print(pose.residue(167).name3())

    def setup_specialrots(self, pose):
        '''Set up the special rotamers'''
        self.specialrot_dict = {}
        self.specialrot_pose = pose.clone()
        for position in self.special_rotamers:
            current_rsd_type_ptr = pose.residue_type_ptr(position)
            new_rsd_type_mutable = rosetta.core.chemical.MutableResidueType(current_rsd_type_ptr)
            new_rsd_type_mutable.add_variant_type(rosetta.core.chemical.SPECIAL_ROT)
            new_rsd_type = rosetta.core.chemical.ResidueType.make(new_rsd_type_mutable)
            rosetta.core.pose.replace_pose_residue_copying_existing_coordinates(self.specialrot_pose,
                                                                                position, new_rsd_type)
            self.specialrot_dict[position] = self.specialrot_pose.residue(position).clone()

    def accept_to_best(self, pose):
        '''Calculate pose score and remember the pose if it is the best score we've seen'''
        current_score = self.sfxn(pose)
        if current_score < self.best_score:
            self.best_score = current_score
            self.best_pose = pose.clone()
            print(f'Accepted new best score: {self.best_score}', flush=True)

        return self.best_pose

    def setup_minmover(self, tolerance):
        '''Setup default minmover'''
        mintype = 'lbfgs_armijo_nonmonotone'
        minimizer = rosetta.protocols.minimization_packing.MinMover()
        minimizer.movemap(self.movemap)
        minimizer.score_function(self.sfxn)
        minimizer.min_type(mintype)
        minimizer.tolerance(tolerance)
        return minimizer

    def apply(self, pose):
        '''
        Initialize packer, edit rotamerset
        Loop through packing & min
        '''

        if not self.taskfactory:
            self.setup_default_taskfactory()

        self.parse_script(self.script)
        if not self.movemap:
            self.movemap = self.setup_default_movemap(pose)

        # Setup special rotamers
        self.setup_specialrots(pose)
        self.best_pose = pose

        # Setup constraints
        if self.bb_cst:
            coord_cst = constraint_generator.CoordinateConstraintGenerator()
            if not self.sc_cst:
                coord_cst.set_sidechain(False)
            constraints = coord_cst.apply(pose)
            for cst in constraints:
                pose.add_constraint(cst)

        # Start of protocol
        for outer in range(0, self.nrepeats):
            for cmd in self.script_steps:
                if cmd[0] == 'coord_cst_weight':
                    if self.ramp_down_constraints:
                        self.sfxn.set_weight(ScoreType.coordinate_constraint, cmd[1])
                elif cmd[0].startswith('scale'):
                    st = cmd[0].split(':')[1]
                    if st == 'fa_rep':
                        self.sfxn.set_weight(ScoreType.fa_rep, cmd[1])
                    if st == 'fa_atr':
                        self.sfxn.set_weight(ScoreType.fa_atr, cmd[1])
                elif cmd[0] == 'repack':
                    self.pack_rotamers(pose)
                elif cmd[0] == 'min':
                    minmover = self.setup_minmover(cmd[1])
                    minmover.apply(pose)
                elif cmd[0] == 'accept_to_best':
                    pose = self.accept_to_best(pose)


def test_parse_script():
    init()
    des = SpecialRotDesign()
    print(des.script)
    des.parse_script(des.script)
    print(des.script_steps)


def test_fastdesign():
    init()
    pdbfile = "/Users/codykrivacic/intelligent_design/helix_workspaces/002_1b8d_4_3_16.pdb.gz"
    ref = create_score_function('ref2015')
    pose = pose_from_file(pdbfile)
    print(f"Initial score: {ref(pose)}")
    des = SpecialRotDesign(special_rotamers=[167], special_rot_weight=-1.5, ramp_down_constraints=True)
    des.apply(pose)
    final = des.best_pose
    ref(final)
    final.dump_pdb('/Users/codykrivacic/intelligent_design/helix_workspaces/testout.pdb.gz')
    print(f"Final score:{ref(final)}")


def test_rotamer_sets():
    pdbfile = "/Users/codykrivacic/intelligent_design/helix_workspaces/002_1b8d_4_3_16.pdb.gz"
    init()
    pose = pose_from_file(pdbfile)
    des = SpecialRotDesign(special_rotamers=[167], special_rot_weight=-1.5)
    des.setup_specialrots(pose)
    des.pack_rotamers(pose)


if __name__ == '__main__':
    # test_parse_script()
    # test_rotamer_sets()
    test_fastdesign()