'''
Design the interface and surroundings of a matched scaffold.
Requires score dataframe and match dataframe.

Usage:
    design_interface.py <workspace> [options]

Options:
    --nstruct=INT, -n  How many designs to make for each input
    --task=INT  Run  a specific task
    --special-rot  Use special rotamer bonus to keep good residues instead of turning off packing
    --suffix=STR  Add a suffix to all output names
'''
import sys, os, json
import pandas as pd
import numpy as np
from helix import workspace as ws
from helix import big_jobs
from helix.utils.numeric import euclidean_distance
from helix.patchman.pythondesign import SpecialRotDesign
from helix.rifdock import interface
from helix.utils import utils
from helix.patchman.design_patchman import select_good_residues
import docopt
# Basic PyRosetta stuff
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
from pyrosetta import rosetta
from pyrosetta.rosetta.core.scoring import CA_rmsd
from pyrosetta.rosetta.core.scoring import all_atom_rmsd
# Packing stuff
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
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


def get_layer_design():
    # Not used

   '''
    <RESIDUE_SELECTORS>
        # These are the default paramters to the old TaskOperation LayerDesign
        <Layer name="surface" select_core="false" select_boundary="false" select_surface="true" use_sidechain_neighbors="true"/>
        <Layer name="boundary" select_core="false" select_boundary="true" select_surface="false" use_sidechain_neighbors="true"/>
        <Layer name="core" select_core="true" select_boundary="false" select_surface="false" use_sidechain_neighbors="true"/>
        <SecondaryStructure name="sheet" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="E"/>
        <SecondaryStructure name="entire_loop" overlap="0" minH="3" minE="2" include_terminal_loops="true" use_dssp="true" ss="L"/>
        <SecondaryStructure name="entire_helix" overlap="0" minH="3" minE="2" include_terminal_loops="false" use_dssp="true" ss="H"/>
        <And name="helix_cap" selectors="entire_loop">
            <PrimarySequenceNeighborhood lower="1" upper="0" selector="entire_helix"/>
        </And>
        <And name="helix_start" selectors="entire_helix">
            <PrimarySequenceNeighborhood lower="0" upper="1" selector="helix_cap"/>
        </And>
        <And name="helix" selectors="entire_helix">
            <Not selector="helix_start"/>
        </And>
        <And name="loop" selectors="entire_loop">
            <Not selector="helix_cap"/>
        </And>

    </RESIDUE_SELECTORS>

    <TASKOPERATIONS>
        <DesignRestrictions name="layer_design_F_boundary_M"> # 2 / 3 -- This is 3 / 3 if you don't have approx, StructProfile, and sap
                                                                         Rosetta will put hydrophobics on your surface and polars in your core
                                                                         We do need core polars for interfaces though
            <Action selector_logic="surface AND helix_start"  aas="DEHKPQR"/>
            <Action selector_logic="surface AND helix"        aas="EHKQR"/>
            <Action selector_logic="surface AND sheet"        aas="EHKNQRST"/>
            <Action selector_logic="surface AND loop"         aas="DEGHKNPQRST"/>
            <Action selector_logic="boundary AND helix_start" aas="ADEFHIKLMNPQRSTVWY"/>
            <Action selector_logic="boundary AND helix"       aas="ADEFHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND sheet"       aas="DEFHIKLMNQRSTVWY"/>
            <Action selector_logic="boundary AND loop"        aas="ADEFGHIKLNPQRSTVWY"/>
            <Action selector_logic="core AND helix_start"     aas="AFILMPVWY"/>
            <Action selector_logic="core AND helix"           aas="AFILMVWYDENQSTH"/>
            <Action selector_logic="core AND sheet"           aas="FILMVWYDENQSTH"/>
            <Action selector_logic="core AND loop"            aas="AFGILPVWYDENQSTH"/>
            <Action selector_logic="helix_cap"                aas="DNST"/>
        </DesignRestrictions>
    </TASKOPERATIONS>
   '''
   layer_design = XmlObjects.static_get_task_operation('layer_design_F_boundary_M')
   return layer_design


def apply_filters(workspace, pose, input_pose=None):
    psipred_single = os.path.expanduser('~/software/fragments/psipred/runpsipred_single')
    filter_objs = {}
    ss_vall = workspace.find_path('ss_grouped_vall_all.h5')
    if os.path.exists(ss_vall):
        worst_9mer_filters = f'''
        <FILTERS>
            <worst9mer name="worst_9mer" confidence="0" rmsd_lookup_threshold="0.01" report_mean_median="true" />
            <worst9mer name="worst_9mer_helix" confidence="0" rmsd_lookup_threshold="0.01" report_mean_median="true" 
            only_helices="true" />
        </FILTERS>
        '''
        # fragments_xml = XmlObjects.create_from_string(worst_9mer_filters)
        # for filter_name in ['worst_9mer', 'worst_9mer_helix']:
        #     filter_objs[filter_name] = fragments_xml.get_filter(filter_name)
    filters_string = f'''
    <FILTERS>
        <ResidueCount name="res_count_all" max_residue_count="9999" confidence="0"/>
        <ResidueCount name="ala_count" residue_types="ALA" max_residue_count="6" confidence="0"/>
        <SSPrediction name="mismatch_probability" confidence="0" cmd="{psipred_single}" use_probability="1" mismatch_probability="1" use_svm="0" />
        <SSShapeComplementarity name="ss_sc" verbose="0" confidence="0" min_sc="0.800" />
        <Time name="time"/>
    </FILTERS>
    '''
    filter_xml = XmlObjects.create_from_string(filters_string)
    for filter_name in ['res_count_all', 'ala_count', 'mismatch_probability',
                        'ss_sc', 'time']:
        filter_objs[filter_name] = filter_xml.get_filter(filter_name)

    sap = '''
    <RESIDUE_SELECTORS>
        <Chain name="chainA" chains="A"/>
        <Chain name="chainB" chains="B"/>
        <True name="true_sel" />
    </RESIDUE_SELECTORS>
    <SIMPLE_METRICS>

        <SapScoreMetric name="sap_score" score_selector="chainA" />
        <SapScoreMetric name="sap_score_target" score_selector="chainB" />
        <SapScoreMetric name="binder_blocked_sap" score_selector="chainA" sap_calculate_selector="chainA" sasa_selector="true_sel" />
        <SapScoreMetric name="target_blocked_sap" score_selector="chainB" sap_calculate_selector="chainB" sasa_selector="true_sel" />

        <CalculatorMetric name="binder_delta_sap" equation="binder_sap_score - binder_blocked_sap" >
            <VAR name="binder_sap_score" metric="sap_score"/>
            <VAR name="binder_blocked_sap" metric="binder_blocked_sap"/>
        </CalculatorMetric>

        <CalculatorMetric name="target_delta_sap" equation="target_sap_score - target_blocked_sap" >
            <VAR name="target_sap_score" metric="sap_score_target"/>
            <VAR name="target_blocked_sap" metric="target_blocked_sap"/>
        </CalculatorMetric>

    </SIMPLE_METRICS>
    '''
    metric_objs = {}
    metrics_xml = XmlObjects.create_from_string(sap)
    metric_objs['sap_score'] = metrics_xml.get_simple_metric('sap_score')
    metric_objs['sap_score_target'] = metrics_xml.get_simple_metric('sap_score_target')
    metric_objs['binder_delta_sap'] = metrics_xml.get_simple_metric('binder_delta_sap')
    metric_objs['target_delta_sap'] = metrics_xml.get_simple_metric('target_delta_sap')

    buns_all = '''
    <BuriedUnsatHbonds 
        name="Buried Unsat [[-]]"
        report_all_heavy_atom_unsats="true" scorefxn="ref2015"
        cutoff="4" residue_surface_cutoff="20.0"
        ignore_surface_res="true" print_out_info_to_pdb="true"
        dalphaball_sasa="1" probe_radius="1.1" confidence="0"
        only_interface="true" />

    '''
    buns_interface = '''
    <RESIDUE_SELECTORS>
        <Chain name="chainA" chains="A"/>
        <Chain name="chainB" chains="B"/>
        <Neighborhood name="interface_chA" selector="chainB" distance="10.0" /> # 1 / 3 -- interface size. 8 - 12 seems reasonable
        <Neighborhood name="interface_chB" selector="chainA" distance="10.0" /> #          at 8, you have trouble with ARG and LYS though
        <And name="AB_interface" selectors="interface_chA,interface_chB" />
    </RESIDUE_SELECTORS>
    <FILTERS>
        <BuriedUnsatHbonds name="buns_interface" use_reporter_behavior="true" report_all_heavy_atom_unsats="true" scorefxn="sfxn" residue_selector="AB_interface" ignore_surface_res="false" print_out_info_to_pdb="true" confidence="0" use_ddG_style="true" burial_cutoff="0.01" dalphaball_sasa="true" probe_radius="1.1" max_hbond_energy="1.5" burial_cutoff_apo="0.2" />
    </FILTERS>
    '''
    buns_sc = '''
    <BuriedUnsatHbonds 
        name="Buried Unsat Sidechains [[-]]"
        report_sc_heavy_atom_unsats="true" scorefxn="ref2015" cutoff="4"
        residue_surface_cutoff="20.0" ignore_surface_res="true"
        print_out_info_to_pdb="true"  dalphaball_sasa="1"
        probe_radius="1.1" confidence="0" only_interface="true" />

    '''
    npsa = '''
    <BuriedSurfaceArea name="Buried Nonpolar Surface Area [[+]]"
    select_only_FAMILYVW="true" filter_out_low="false"
    atom_mode="all_atoms"
    confidence="1.0"
    />
    '''
    exposed_hydrophobics = \
        '''
    <ExposedHydrophobics
    name="ExposedHydrophobics SASA [[-]]"
    sasa_cutoff="20"
    threshold="-1"
    />
    '''
    packstat = '''
      <PackStat
        name="PackStat Score [[+]]"
        threshold="0"
      />
    '''
    sc = \
        '''
    <ShapeComplementarity name="shapecomp" min_sc="0" 
    jump="1"  write_int_area="false" />
    '''
    ia = \
        '''
    <InterfaceAnalyzerMover name="interface_analyzer" 
    pack_separated="false" pack_input="false"
    resfile="false" packstat="true"
    interface_sc="false" tracer="false"
    use_jobname="false" 
    interface="A_B" />
    '''
    contact = \
        '''
    <RESIDUE_SELECTORS>
        <Chain name="chA" chains="A"/>
        <Chain name="chB" chains="B"/>
    </RESIDUE_SELECTORS>
    <FILTERS>
        <ContactMolecularSurface name="contact" target_selector="chA"
        binder_selector="chB" />
    </FILTERS>
    '''
    ia_mover = XmlObjects.static_get_mover(ia)
    ia_mover.set_compute_interface_sc(False)
    # ia_mover.set_compute_separated_sasa(False)
    # ia_mover.set_calc_dSASA(False)

    poseA = utils.pose_get_chain(pose, 'A')
    poseB = utils.pose_get_chain(pose, 'B')
    ref = create_score_function('ref2015')
    ref(poseA)
    ref(poseB)

    filter_objs['buns_all'] = XmlObjects.static_get_filter(buns_all)
    filter_objs['buns_sc'] = XmlObjects.static_get_filter(buns_sc)
    buns_int_xml = XmlObjects.create_from_string(buns_interface)
    filter_objs['buns_interface'] = buns_int_xml.get_filter('buns_interface')
    filter_objs['npsa'] = XmlObjects.static_get_filter(npsa)
    filter_objs['exposed_hydrophobics'] = XmlObjects.static_get_filter(exposed_hydrophobics)
    filter_objs['packstat']= XmlObjects.static_get_filter(packstat)
    filter_objs['shape_complementarity'] = XmlObjects.static_get_filter(sc)
    filter_objs['contact_molecular_surface'] = XmlObjects.create_from_string(contact).get_filter('contact')

    row = {}
    for filter_name in filter_objs:
        row[filter_name] = filter_objs[filter_name].report_sm(pose)
    for metric_name in metric_objs:
        row[metric_name] = metric_objs[metric_name].calculate(pose)

    interface_scorer = interface.InterfaceScore(pose)
    interface_score = interface_scorer.apply()
    row['interface_score'] = interface_score
    row['n_hbonds'] = interface_scorer.n_hbonds

    npsa_complex = row['npsa']
    npsa_A = filter_objs['npsa'].report_sm(poseA)
    npsa_B = filter_objs['npsa'].report_sm(poseB)
    row['delta_npsa'] = npsa_complex - npsa_B - npsa_A

    ia_mover.apply(pose)
    # Can't pickle a C++ set, so put it in a Python list
    int_set = []
    for interface_resi in ia_mover.get_interface_set():
        int_set.append(interface_resi)

    row['n_interface_residues'] = ia_mover.get_num_interface_residues()
    row['complex_sasa'] = ia_mover.get_complexed_sasa()
    row['delta_sasa'] = ia_mover.get_interface_delta_sasa()
    row['crossterm_energy'] = ia_mover.get_crossterm_interface_energy()
    row['interface_packstat'] = ia_mover.get_interface_packstat()
    row['delta_unsat'] = ia_mover.get_interface_delta_hbond_unsat()
    row['interface_dG'] = ia_mover.get_interface_dG()
    row['interfac_residues'] = int_set

    if input_pose:
        ca_rmsd = CA_rmsd(pose, input_pose)
        aa_rmsd = all_atom_rmsd(pose, input_pose)
        row['ca_rmsd'] = ca_rmsd
        row['all_atom_rmsd'] = aa_rmsd

    return row


def calculate_fsf(workspace, pose, insertion, suffix, test_run=False):
    fsf_worst = '''
      <FragmentScoreFilter
        name="Worst. 9-Residue Fragment Crmsd ({largest_loop_start})"
        scoretype="FragmentCrmsd"
        sort_by="FragmentCrmsd"
        threshold="9999" 
        direction="-"
        start_res="{largest_loop_start}"
        end_res="{largest_loop_end}"
        compute="maximum"
        outputs_folder="{seqprof_dir}"
        outputs_name="{suffix}" 
        csblast="/wynton/home/kortemme/krivacic/software/fragments/csblast-2.2.3_linux64"  
        blast_pgp="/wynton/home/kortemme/krivacic/software/fragments/blast/blast-2.2.26/bin/blastpgp" 
        psipred="/wynton/home/kortemme/krivacic/software/fragments/psipred/runpsipred_single" 
        sparks-x="/wynton/home/kortemme/krivacic/software/fragments/sparks-x" 
        sparks-x_query="/wynton/home/kortemme/krivacic/software/fragments/sparks-x/bin/buildinp_query.sh" 
        frags_scoring_config="{fragment_weights_path}"
        placeholder_seqs="/wynton/home/database/blast/blastdb/pdbaa"
        print_to_pdb="true"
        n_frags="200"
        n_candidates="1000" 
        fragment_size="9"
        vall_path="{vall_path}"
      />
    '''.format(largest_loop_start=insertion['start'],
               largest_loop_end=insertion['stop'],
               seqprof_dir=workspace.seqprof_dir, suffix=suffix,
               fragment_weights_path=workspace.fragment_weights_path,
               vall_path=workspace.rosetta_vall_path(test_run))

    fsf_obj = XmlObjects.static_get_filter(fsf_worst)
    try:
        # First try to run FSF because it needs to generate the
        # other sequence profile files
        score = fsf_obj.report_sm(pose)
    except:
        import subprocess
        # Now get rid of the empty .fasta.pssm file it creates
        # and try again in Python
        rempath = os.path.join(workspace.seqprof_dir,
                               '{}.fasta.phipsi'.format(suffix))
        print('REMOVING {}'.format(rempath))
        os.remove(rempath)
        print('CWD: {}'.format(os.getcwd()))
        print('TASK: {}'.format(suffix))
        print('FILES IN CWD:')
        print(os.listdir(os.getcwd()))
        cmd = [
            '/wynton/home/kortemme/krivacic/software/fragments/sparks-x/bin/buildinp_query.sh',
            os.path.join(workspace.seqprof_dir, '{}.fasta'.format(suffix)),
        ]

        process = subprocess.run(cmd,
                                 env=dict(SPARKSXDIR='/wynton/home/kortemme/krivacic/software/fragments/sparks-x',
                                          **os.environ))
        score = fsf_obj.report_sm(pose)

    return score


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
    '''
    Class to help manage things for interface design.
    If special_rot is set to True, special rotamer bonuses will be applied as normal. Otherwise, good rotamers
    will simply be frozen. Rotamers will only be considered "good" if they are comparable to natural
    rotamers in similar environments in terms of cross-chain score AND were transferred from a docked helix.
    '''
    def __init__(self, workspace, df_group, task_id, special_rot=False, special_rot_weight=-3.0, buns_penalty=True,
                 ramp_cst=True, test_run=False, interface_upweight=True, total_inputs=0,
                 suffix=False):
        self.workspace = workspace

        self.df = df_group
        self.task_id = task_id
        self.pdb_path = os.path.join(self.workspace.root_dir, self.df.iloc[0].superimposed_file)

        # Initiate Rosetta
        dalphaball = os.path.join(self.workspace.rosetta_dir,
                                  'source', 'external', 'DAlpahBall',
                                  'DAlphaBall.gcc')
        ss_vall = self.workspace.find_path('ss_grouped_vall_all.h5')
        init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
             ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false ' \
             '-indexed_structure_store:fragment_store {} ' \
             '-in:file:s {}'.format(dalphaball, ss_vall, self.pdb_path))

        self.design_pose = pose_from_file(
            self.pdb_path
        )
        self.input_pose = self.design_pose.clone()
        sfxn = create_score_function('ref2015')
        sfxn.set_weight(ScoreType.coordinate_constraint, 1.0)
        sfxn.set_weight(ScoreType.aa_composition, 1.0)
        sfxn.set_weight(ScoreType.arg_cation_pi, 3.0)

        self.summarized_residue_scores = utils.safe_load(workspace.find_path(
            'summarized_res_scores.pkl'
        ))
        self.special_rot = special_rot

        if buns_penalty:
            # If buried unsat penalty, use this sfxn
            buns_sfxn = '''
            <ScoreFunction name="sfxn" weights="ref2015" >
                <Reweight scoretype="approximate_buried_unsat_penalty" weight="5.0" />
                <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.5" />
                <Set approximate_buried_unsat_penalty_burial_atomic_depth="3.5" />
                <Set approximate_buried_unsat_penalty_hbond_bonus_cross_chain="-1" />
            </ScoreFunction>
            '''
            sfxn = XmlObjects.static_get_score_function(buns_sfxn)
        if special_rot:
            sfxn.set_weight(rosetta.core.scoring.special_rot,
                            special_rot_weight)
        self.sfxn_cst = sfxn

        self.ramp_cst = ramp_cst
        self.interface_upweight=interface_upweight
        self.test_run=test_run
        basename = os.path.basename(self.pdb_path).split('.')[0] + f'_{self.task_id//total_inputs}'
        if suffix:
            basename += suffix
            self.suffix = suffix
        else:
            self.suffix = ''
        if not os.path.exists(self.workspace.design_dir):
            os.makedirs(self.workspace.design_dir, exist_ok=True)
        self.output_pickle = os.path.join(self.workspace.design_dir, basename + '.pkl')
        basename += '.pdb.gz'
        self.output_file = os.path.join(self.workspace.design_dir, basename)


    def get_json(self):
        model_no = os.path.basename(self.pdb_path).split('_')[1]
        lhl_folder = os.path.join(self.workspace.root_dir, '..', 'regenerated_data_sets_2020_03',
                                  'sequence_design_for_LHL_reshaping_2lv8_two_LHL',
                                  'selected_designs_for_state_count')
        insertion_file = os.path.join(lhl_folder, f'insertion_points_{model_no}.json')
        with open(insertion_file, 'r') as f:
            insertions = json.load(f)
        return insertions

    def filter(self):
        # row = {}
        row = apply_filters(self.workspace, self.design_pose, self.input_pose)
        row['superimposed_file'] = self.df.iloc[0]['superimposed_file']
        row['design_file'] = os.path.relpath(self.output_file, self.workspace.root_dir)
        row['suffix'] = self.suffix
        ref = create_score_function('ref2015')
        row['total_score'] = ref(self.design_pose)
        # self.design_pose.dump_pdb(self.output_file)
        # self.row.to_pickle(self.output_pickle)

        i = 0
        for insertion in self.get_json():
            i += 1
            try:
                row[f"frag_score_filter_{i}"] = calculate_fsf(self.workspace, self.design_pose, insertion,
                                                            f"{self.suffix}_{str(self.task_id)}_{i}",
                                                            test_run=self.test_run)
                                                            # test_run=True)
            except:
                row[f"frag_score_filter_{i}"] = np.nan
        self.row = row

    def setup_relax_task_factory(self):
        tf_str = '''
        <RESIDUE_SELECTORS>
            <InterfaceByVector name="interface_selector" cb_dist_cut="12" nearby_atom_cut="6.5" vector_angle_cut="75">
               <Chain chains='A'/>
               <Chain chains='B'/>
            </InterfaceByVector>
            <Chain name="chainB" chains="B"/>
            <Not name="not_interface"  selector="interface_selector" />
            <And name="chainB_not_interface" selectors="not_interface,chainB" />
        </RESIDUE_SELECTORS>
        <TASKOPERATIONS>
            <OperateOnResidueSubset name="restrict_target_not_interface" selector="chainB_not_interface">
                <PreventRepackingRLT/>
            </OperateOnResidueSubset>
            <ExtraRotamersGeneric name="ex1_ex2" ex1="1" ex2aro="1" ex2="0" />
            <LimitAromaChi2 name="limitchi2" chi2max="110" chi2min="70" include_trp="True" />
            <PruneBuriedUnsats name="prune_buried_unsats" allow_even_trades="false" atomic_depth_cutoff="3.5" minimum_hbond_energy="-0.5" />
        </TASKOPERATIONS>
        '''
        tf_obj = XmlObjects.create_from_string(tf_str)
        tf = TaskFactory()
        for taskop in ['restrict_target_not_interface', 'ex1_ex2', 'limitchi2', 'prune_buried_unsats']:
            tf.push_back(tf_obj.get_task_operation(taskop))

        return tf


    def setup_design_task_factory(self, initial_design=False):
        '''Set up task factory for design'''
        selector = residue_selector.ChainSelector('A')

        interface_selector_str = \
        '''
        <RESIDUE_SELECTORS>
            <InterfaceByVector name="interface_vector" cb_dist_cut="11" nearby_atom_cut="6.5" vector_angle_cut="75">
               <Chain chains='A'/>
               <Chain chains='B'/>
            </InterfaceByVector>
            <Chain name="chainA" chains="A"/>
            <Chain name="chainB" chains="B"/>
            <Neighborhood name="interface_chA" selector="chainB" distance="10.0" /> # 1 / 3 -- interface size. 8 - 12 seems reasonable
            <Neighborhood name="interface_chB" selector="chainA" distance="10.0" /> #          at 8, you have trouble with ARG and LYS though
            <And name="AB_interface" selectors="interface_chA,interface_chB" />
            <Or name="any_interface" selectors="AB_interface,interface_vector" />
            <Not name="Not_interface" selector="AB_interface" />
            <And name="actual_interface_chA" selectors="AB_interface,chainA" />
            <And name="actual_interface_chB" selectors="AB_interface,chainB" />
            <And name="chainB_not_interface" selectors="Not_interface,chainB" />
        </RESIDUE_SELECTORS>
        '''
        interface_selector_xml = XmlObjects.create_from_string(interface_selector_str)
        interface_selector = interface_selector_xml.get_residue_selector('interface_vector')

        include_current = XmlObjects.static_get_task_operation(
            '''<IncludeCurrent name="include_current" />'''
        )
        limit_chi2 = XmlObjects.static_get_task_operation(
            '''<LimitAromaChi2 name="limit_chi2" chi2max="110" chi2min="70" include_trp="True" />'''
        )
        exchi = XmlObjects.static_get_task_operation(
            '''<ExtraRotamersGeneric name="extra_rot" ex1="1" ex2aro="1" ex2="0" />'''
        )
        no_gly = XmlObjects.static_get_task_operation(
            '''<DisallowIfNonnative name="no_gly" resnum="0" disallow_aas="G" />'''
        )
        no_pro = XmlObjects.static_get_task_operation(
            '''<DisallowIfNonnative name="no_pro" resnum="0" disallow_aas="P" />'''
        )

        tf = TaskFactory()
        tf.push_back(operation.InitializeFromCommandline())
        tf.push_back(include_current)
        tf.push_back(limit_chi2)
        tf.push_back(exchi)
        tf.push_back(no_gly)
        tf.push_back(no_pro)

        no_packing = operation.PreventRepackingRLT()
        no_design = operation.RestrictToRepackingRLT()

        if initial_design:
            nopack = self.special_residues
        else:
            self.get_good_residues()
            nopack = self.nopack
            print('GOOD RESIDUES: ')
            print(nopack)
        if not self.special_rot and not initial_design:
            # If not using special rot score term, freeze the "good rotamers" in place
            if len(self.nopack) > 0:
                nopack_selector = utils.list_to_res_selector(nopack)
                good_rots = operation.OperateOnResidueSubset(no_packing,
                                                             nopack_selector)
                tf.push_back(good_rots)
            else:
                nopack_selector = residue_selector.FalseResidueSelector()
        elif initial_design:
            nopack_selector = utils.list_to_res_selector(nopack)
            good_rots = operation.OperateOnResidueSubset(no_design, nopack_selector)
            tf.push_back(good_rots)
        else:
            nopack_selector = residue_selector.FalseResidueSelector()

        # or_selector = residue_selector.OrResidueSelector(selector,
        #                                                  interface_selector)

        if initial_design:
            clash_selector_input = residue_selector.OrResidueSelector(interface_selector, nopack_selector)
        else:
            clash_selector_input = interface_selector
        clash_selector = rosetta.core.pack.task.residue_selector.ClashBasedShellSelector(clash_selector_input)
        clash_selector.set_num_shells(2)
        clash_selector.invert(False)
        clash_selector.set_include_focus(True)
        interface_residues = utils.res_selector_to_size_list(interface_selector.apply(self.design_pose))
        clash_resis = utils.res_selector_to_size_list(clash_selector.apply(self.design_pose), pylist=True)
        print('INTERFACE SELECTOR:')
        print(' or resi '.join([str(x) for x in interface_residues]))
        print('CLASH SELECTOR:')
        print(' or resi '.join([str(x) for x in clash_resis]))

        if self.interface_upweight:
            upweight = \
            '''
            <ProteinProteinInterfaceUpweighter name="upweight_interface" interface_weight="2" />
            '''
            upweight_taskop = XmlObjects.static_get_task_operation(upweight)
            tf.push_back(upweight_taskop)

        # Interface and chain B
        and_selector = residue_selector.AndResidueSelector(interface_selector,
                                                           selector)
        not_interface = residue_selector.NotResidueSelector(and_selector)
        notdesign = operation.OperateOnResidueSubset(no_design,
                                                     not_interface)

        # Clash-based repack shell (including focus) and nopack selector
        or_selector = residue_selector.OrResidueSelector(clash_selector, nopack_selector)
        not_selector = residue_selector.NotResidueSelector(or_selector)
        static = operation.OperateOnResidueSubset(no_packing,
                                                  not_selector)

        # tf.push_back(notaa)
        tf.push_back(notdesign)
        tf.push_back(static)

        return tf

    def design(self):
        '''Design interface residues after transferring from helix'''
        # Set to interface design script
        interface_script_path = os.path.join(
            self.workspace.rosetta_dir, 'database', 'sampling',
            'relax_scripts', 'InterfaceDesign2019.txt'
        )
        sap_cst_str = '''
        <RESIDUE_SELECTORS>
            <Chain name="chainA" chains="A"/>
        </RESIDUE_SELECTORS>
        <MOVERS>
            <AddSapConstraintMover name="add_sap" speed="lightning" sap_goal="0" penalty_per_sap="1" score_selector="chainA" sap_calculate_selector="chainA" /> 
            <AddCompositionConstraintMover name="3trp" > # 2 / 3 -- penalize TRP by 3. This leads to a reasonable number of TRP given everything else
                <Comp entry="PENALTY_DEFINITION;TYPE TRP;ABSOLUTE 0;PENALTIES 0 3;DELTA_START 0;DELTA_END 1;BEFORE_FUNCTION CONSTANT;AFTER_FUNCTION LINEAR;END_PENALTY_DEFINITION;" />
            </AddCompositionConstraintMover>
            <AddCompositionConstraintMover name="2met" > # 2 / 3 -- Rosetta loves MET, but the Bakerlab doesn't.
                <Comp entry="PENALTY_DEFINITION;TYPE MET;ABSOLUTE 0;PENALTIES 0 2;DELTA_START 0;DELTA_END 1;BEFORE_FUNCTION CONSTANT;AFTER_FUNCTION LINEAR;END_PENALTY_DEFINITION;" />
            </AddCompositionConstraintMover>
        </MOVERS>
        '''
        sap_cst_xml = XmlObjects.create_from_string(sap_cst_str)
        sap_cst_mvr = sap_cst_xml.get_mover("add_sap")
        sap_cst_mvr.apply(self.design_pose)
        trp_mvr = sap_cst_xml.get_mover('3trp')
        trp_mvr.apply(self.design_pose)
        met_mvr = sap_cst_xml.get_mover('2met')
        met_mvr.apply(self.design_pose)
        movemap = self.setup_design_movemap()

        tf_initial = self.setup_design_task_factory(initial_design=True)
        # self.design_pose.dump_pdb('testout_selectors.pdb')
        # Quick FastDesign ( 1 repeat ) just to get things in the right place so we can determine if we want
        # to keep special rotamers.
        # fastdes_initial = rosetta.protocols.denovo_design.movers.FastDesign(self.sfxn_cst,
        #                                                                     interface_script_path)
        fastdes_initial = rosetta.protocols.denovo_design.movers.FastDesign(self.sfxn_cst,
                                                                    1)
        fastdes_initial.set_task_factory(tf_initial)
        fastdes_initial.set_movemap(movemap)
        print('Performing initial design')
        fastdes_initial.apply(self.design_pose)
        # self.design_pose.dump_pdb('test_design.pdb')

        tf_final = self.setup_design_task_factory(initial_design=False)
        # Get rid of constraints and only add back those that are determined to be "good residues" after initial fastdesign
        self.design_pose.remove_constraints()
        for resi in self.nopack:
            print(f'Adding sidechain constraints for residue {resi}')
            for cst in self.sc_constraints[resi]:
                self.design_pose.add_constraint(cst)
        self.design_pose.add_constraint(self.bb_constraints)


        if not self.special_rot:
            fastdes = rosetta.protocols.denovo_design.movers.FastDesign(self.sfxn_cst,
                                                                                  interface_script_path)
        else:
            fastdes = SpecialRotDesign(special_rotamers=self.nopack,
                                       bb_cst=True, rosettadir=self.workspace.rosetta_dir,
                                       script=interface_script_path, taskfactory=tf_final)
            fastdes.set_scorefxn(self.sfxn_cst)
            fastdes.set_movemap(movemap)
        if self.ramp_cst:
            fastdes.ramp_down_constraints(True)
        else:
            fastdes.ramp_down_constraints(False)

        fastdes.set_task_factory(tf_final)
        for i in range(0, 2):
            print(f'Performing FastDesign round {i}')
            fastdes.apply(self.design_pose)

        # Relax w/o constraints
        self.design_pose.remove_constraints()
        self.design_pose.clear_sequence_constraints()

        fastrelax_str = '''
        <MOVERS>
            <FastRelax name="FastRelax" repeats="1" batch="false" ramp_down_constraints="false" 
            cartesian="false" bondangle="false" bondlength="false" min_type="dfpmin_armijo_nonmonotone" >
                <MoveMap name="MM" >
                    <Chain number="1" chi="true" bb="true" />
                    <Chain number="2" chi="true" bb="false" />
                    <Jump number="1" setting="true" />
                </MoveMap>
            </FastRelax>
        </MOVERS>
        '''
        fastrelax_xml = XmlObjects.create_from_string(fastrelax_str)
        fastrelax_mover = fastrelax_xml.get_mover('FastRelax')
        fastrelax_mover.set_task_factory(self.setup_relax_task_factory())
        fastrelax_mover.set_movemap(movemap)
        fastrelax_mover.set_scorefxn(self.sfxn_cst)
        fastrelax_mover.apply(self.design_pose)

    def get_good_residues(self):
        '''Find good rotamers as compared to natural protein interfaces'''
        # No need to save the pose because we already have a minimized pose
        self.nopack, pose = select_good_residues(self.design_pose, self.summarized_residue_scores, is_pose=True,
                                                          cst_sc=True, minimize=False, chain='A')
        print('ORIGINAL NOPACK RESIDUES: ', self.nopack)
        print('SPECIAL RESIDUES: ', self.special_residues)
        self.nopack = [x for x in self.nopack if x in self.special_residues]
        print("NOPACK RESIDUES: ", self.nopack)

    def transfer_residue(self, pose2, pose1_resnum, pose2_resnum):
        '''Transfers a rotamer from pose2 to pose1'''
        print(f"Transferring residue {pose2_resnum}, {pose2.residue(pose2_resnum).name3()} from helix to design position {pose1_resnum}.")
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

        # Now add and save constraints
        self.sc_constraints[pose1_resnum] = []
        for atom in range(pose2_residue.first_sidechain_atom(), pose2_residue.natoms()):
            if pose2_residue.atom_is_hydrogen(atom):
                continue
            xyz = pose2_residue.xyz(atom)
            id = AtomID(atom, int(pose1_resnum))
            reference = AtomID(1, 1)
            func = HarmonicFunc(0, 1)
            cst = CoordinateConstraint(id, reference, xyz, func)
            self.sc_constraints[pose1_resnum].append(cst)
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

        # For design, we'll use neighbor selectors for defining the interface to make sure we capture everything.
        # Here, we want to make sure to only transfer over residues that are pointing towards the interface,
        # otherwise we might transfer a residue to the back of a helix, and during the initial "equilibration" fastdesign
        # (where surrounding residues are allowed to mutate to accommodate the residues from the helices), we might
        # accidentally mess up the protein.
        interface_selector_str =\
        '''
        <InterfaceByVector name="interface_selector" cb_dist_cut="11" nearby_atom_cut="5.5" vector_angle_cut="75">
           <Chain chains='A'/>
           <Chain chains='B'/>
        </InterfaceByVector>
        '''
        interface_selector = XmlObjects.static_get_residue_selector(interface_selector_str)
        selector = residue_selector.ChainSelector('A')
        and_selector = residue_selector.AndResidueSelector(selector,
                                                         interface_selector)

        interface_residues = and_selector.apply(self.design_pose)
        interface_residues = utils.res_selector_to_size_list(interface_residues, pylist=True)
        final_repack_resis = []
        self.sc_constraints = {}
        for idx, row in self.df.iterrows():
            # Helix residue positions are for just that chain
            docked_pose = pose_from_file(
                os.path.join(
                    self.workspace.root_dir, row.design_file
                )
            )
            helix_pose = docked_pose.split_by_chain(2)

            print('')
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
                    if not closest_row.res2 in interface_residues:
                        continue
                    # print(closest_row.res2)
                    # print(closest_row.dist)
                    if closest_row.dist < 1.5:
                        # Relatively relaxed threshold because it might minimize into position
                        self.transfer_residue(helix_pose, closest_row.res2, helix_index)
                        final_repack_resis.append(int(closest_row.res2))


        # Now that all residues have been transferred, repack and minimize
        # For task factory, repack all the interface residues and a clash-based repack shell. Use constraints.
        self.special_residues = final_repack_resis
        print('Constraints present for the following residues:')
        print(', '.join([str(x) for x in final_repack_resis]))
        if len(final_repack_resis) > 0:
            idx_selector = utils.list_to_res_selector(final_repack_resis)
            clash_selector = ClashBasedShellSelector()
            clash_selector.set_focus(idx_selector)
            clash_selector.set_include_focus(True)
            clash_selector.set_num_shells(2)
            clash_selector.invert(True)


            tf = TaskFactory()
            true_selector = residue_selector.TrueResidueSelector()
            no_design = operation.RestrictToRepackingRLT()
            no_des = operation.OperateOnResidueSubset(no_design, true_selector)
            no_packing = operation.PreventRepackingRLT()
            static = operation.OperateOnResidueSubset(no_packing,
                                                      clash_selector)
            tf.push_back(no_des)
            tf.push_back(static)

            # Repack with constraints
            packertask = tf.create_task_and_apply_taskoperations(self.design_pose)
            mover = PackRotamersMover(self.sfxn_cst, packertask)
            mover.apply(self.design_pose)

        # Add backbone constraints for minimization
        coord_cst = constraint_generator.CoordinateConstraintGenerator()
        coord_cst.set_sidechain(False)
        # if not self.sc_cst:
        #     coord_cst.set_sidechain(False)
        self.bb_constraints = coord_cst.apply(self.design_pose)
        for cst in self.bb_constraints:
            self.design_pose.add_constraint(cst)

        packertask = tf.create_task_and_apply_taskoperations(self.design_pose)
        mover = PackRotamersMover(self.sfxn_cst, packertask)
        mover.apply(self.design_pose)

        # Minimize with constraints
        # movemap = self.setup_design_movemap()
        # minmover = MinMover()
        # minmover.movemap(movemap)
        # minmover.score_function(self.sfxn_cst)
        # minmover.apply(self.design_pose)


        # self.design_pose.dump_pdb('testout.pdb')

    def setup_movemap(self):
        # Deprecated
        movemap = MoveMap()
        movemap.set_bb(False)
        movemap.set_chi(False)
        for idx, row in self.df.iterrows():
            start = row.rosetta_helix_start
            stop = row.rosetta_helix_stop
            movemap.set_bb_true_range(start, stop)
            movemap.set_chi_true_range(start, stop)

        return movemap

    def setup_design_movemap(self):
        movemap = MoveMap()
        movemap.set_bb(False)
        movemap.set_chi(False)
        chA_selector = residue_selector.ChainSelector('A')
        chB_selector = residue_selector.ChainSelector('B')
        chA_selection = chA_selector.apply(self.design_pose)
        chB_selection = chB_selector.apply(self.design_pose)
        movemap.set_chi(chA_selection)
        movemap.set_bb(chA_selection)
        movemap.set_chi(chB_selection)
        movemap.set_jump(True)
        return movemap

    def apply(self):
        self.prep_design()
        self.design()
        self.filter()
        self.design_pose.dump_pdb(self.output_file)
        df = pd.DataFrame([self.row])
        df.to_pickle(self.output_pickle)


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
    input_df_path = os.path.join(workspace.complex_dir, 'exported_pdbs.pkl')
    input_df = utils.safe_load(input_df_path)
    inputs = sorted(input_df['superimposed_file'].unique())
    if not job_info['task_id']:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0
    else:
        task_id = int(job_info['task_id'])

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

    designer = InterfaceDesign(workspace, group, task_id, total_inputs=len(inputs), special_rot=args['--special-rot'], suffix=args['--suffix'])
    designer.apply()

if __name__=='__main__':
    # test_prep()
    main()