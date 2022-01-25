"""
Script to run FastDesign + FlexPepDock on outputs of PatchMAN.

Usage:
    design_patchman.py <workspace> [options]

Options:
    --task=INT  Only run a specific task
    --delete  Delete target structures
    --designs-per-task=INT  How many designs per task  [default: 20]
    --align-thresh=FLOAT  Sequence identity score above 
        which favor native residue task operation will be added  [default: 70]
    --buns-penalty  Include a penalty for buried unsat hbonds
    --keep-good-rotamers  Run through positions on docked helix and if
        it is better than the average crosschains core for that residue in
        that environment, keep that rotamer.

"""

from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.protocols import constraint_generator
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
import pyrosetta
from distutils.dir_util import copy_tree
import docopt
import os, sys, glob
import pandas as pd
import pickle5 as pickle
import gzip
from helix import workspace as ws
from helix.utils import utils
from helix import big_jobs
from helix.rifdock import interface
# from helix.matching.scan_helices import contiguous_secstruct
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
import pymol


def strlist_to_vector1_str(strlist):
    # Should also go in utils eventually
    from pyrosetta.rosetta.utility import vector1_std_string
    vector = vector1_std_string()
    for string in strlist:
        vector.append(string)
    return vector


def get_alignment_info(df, path):
    fname = os.path.basename(path)
    row = df[df['complex_basename'] == fname.strip('.gz')]
    length = len(row.iloc[0]['patch_sequence'])
    score = row.iloc[0]['alignment_score']
    patch_seq = row.iloc[0]['patch_sequence'] 
    match_seq = row.iloc[0]['match_sequence']
    matches = 0
    for idx, res in enumerate(patch_seq):
        if res == match_seq[idx]:
            matches += 1
    percent_identity = 100 * (matches / length)
    return score, length, percent_identity


def count_buried_unsat(pose, resnums):
    '''
    See if there are any buried unsatisfied hbonds in a subset of
    residues
    '''
    # sele = utils.list_to_res_selector(resnums)
    buns_all = '''
    <RESIDUE_SELECTORS>
        <Index name="indices" resnums="{}" />
    </RESIDUE_SELECTORS>
    <FILTERS>
    <BuriedUnsatHbonds 
        name="buns_few"
        report_all_heavy_atom_unsats="true" scorefxn="ref2015"
        cutoff="4" residue_surface_cutoff="20.0"
        ignore_surface_res="true" print_out_info_to_pdb="true"
        dalphaball_sasa="1" probe_radius="1.1" confidence="0" 
        residue_selector="indices" />
    </FILTERS>
    '''.format(','.join(resnums))
    # buns_all_obj = XmlObjects.static_get_filter(buns_all)
    buns_all_obj = XmlObjects.create_from_string(buns_all).get_filter('buns_few')
    # buns_all_obj.set_residue_selector(sele)

    return buns_all_obj.report_sm(pose)


def select_good_residues(pdbpath, score_df):
    '''
    Select residues which have similar crosschain scores to natural
    interfaces and which don't have or interact with any buried
    unsatisfied hydrogen bonds. Return a list these resnums,
    with the intent of turning packing off for these residues.
    '''
    from helix.benchmark import score_pdb
    interface = score_pdb.PDBInterface(pdbpath, minimize=True, cst=True)
    interface_scores = interface.interface_all_chains()
    nopack = []
    if interface_scores.shape[0] > 0:
        interface_scores = interface_scores[interface_scores['chain'] == 'B']
        for idx, row in interface_scores.iterrows():
            restype = row['restype']
            burial = row['burial']
            secstruct = row['secstruct']
            score_row = score_df[
                    (score_df['restype']==restype) &
                    (score_df['burial']==burial) &
                    (score_df['secstruct']==secstruct) &
                    (score_df['scoretype']=='total_crosschain')
                    ]
            difference = row['total_crosschain'] - score_row.iloc[0]['median']
            resnum = row['resnum']
            resnums = [str(resnum)]
            for contact in interface.contacts[resnum]:
                resnums.append(str(contact))
            if difference < 0:
                if count_buried_unsat(interface.pose, resnums) < 1:
                    nopack.append(row['resnum'])
                    print("PASSED ALL FILTERS:")
                    print(row)
    else:
        print('No interface found for {}'.format(pdbpath))

    return nopack


def main():
    args = docopt.docopt(__doc__)
    align_threshold = int(args['--align-thresh'])
    try:
        workspace, job_info = big_jobs.initiate()
        # job_info['task_id'] = int(args['--task'])
        # job_info['inputs'] = sorted(glob.glob(
            # os.path.join(workspace.focus_dir, 'patch_*',
                # workspace.scaffold_prefix + '*', 'docked_full',
                # '*.pdb.gz'),
            # ))
        job_info['inputs'] = sorted(glob.glob(
            os.path.join(workspace.focus_dir, 'patch_*',
                workspace.scaffold_prefix + '*', 'docked_full',
                '*.pdb.gz')
            ))
    except:
        print('Maybe this is local?')
        workspace = ws.workspace_from_dir(args['<workspace>'])
        job_info = {
                'task_id': args['--task'],
                'inputs' : sorted(glob.glob(
                    os.path.join(workspace.focus_dir, 'patch_*',
                        workspace.scaffold_prefix + '*', 'docked_full',
                        '*.pdb.gz'),
            ))
                }
    dalphaball = os.path.join(workspace.rosetta_dir,
            'source', 'external', 'DAlpahBall',
            'DAlphaBall.gcc')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro'\
            ' -holes:dalphaball {}'.format(dalphaball))

    if not hasattr(workspace, 'docking_directory'):
        raise Exception("Error: design_patchman.py requires RIFWorkspaces as input. You "\
                "may have provided the root directory to this script "\
                "somehow.")

    inputs = job_info['inputs']
    nstruct = int(args['--designs-per-task'])
    total_jobs = len(inputs)
    print('TOTAL JOBS: {}'.format(total_jobs))

    # print('Job info')
    # print(job_info)
    if not job_info['task_id']:
        if args['--task']:
            task_id = int(args['--task']) - 1
        else:
            task_id = 0
    else:
        task_id = int(job_info['task_id'])
    start = task_id * nstruct
    stop = task_id * nstruct + nstruct
    print(start)
    print(inputs[start])

    target = pymol.cmd.load(workspace.target_path_clean)

    summarized_residue_scores = utils.safe_load(workspace.find_path(
        'summarized_res_scores.pkl'
        ))

    print('TASK: {}'.format(task_id))
    rowlist = []
    latest_patchdir = os.path.dirname(inputs[start])
    print('Starting in directory {}'.format(latest_patchdir))

    alignment_path = os.path.join(latest_patchdir, 'alignment_scores.pkl')
    alignment_df = utils.safe_load(alignment_path)
    alignment_df['complex_basename'] = alignment_df.apply(lambda x:\
            os.path.basename(x['complex']), axis=1)
    for input_idx in range(start, stop):
        # Track patchdir so not always loading alignment score df
        current_patchdir = os.path.dirname(inputs[input_idx])
        if current_patchdir != latest_patchdir:
            alignment_path = os.path.join(current_patchdir,
                    'alignment_scores.pkl')
            alignment_df = utils.safe_load(alignment_path)
            latest_patchdir = current_patchdir
            alignment_df['complex_basename'] = alignment_df.apply(lambda x:\
                    os.path.basename(x['complex']), axis=1)

        # Figure out relative input path for dataframe
        pdb_save = os.path.relpath(inputs[input_idx],
                start=workspace.root_dir)
        pdb = os.path.abspath(inputs[input_idx])
        designed = False

        # Look for residues that have good scores and no BUNS
        nopack = select_good_residues(pdb, summarized_residue_scores)

        # Check if PDB has a score, if so it has been designed already
        # and we can skip
        with gzip.open(pdb, 'rt') as f:
            for line in f:
                if line.startswith('pose'):
                    designed = True

        # Move into pdb folder as working directory
        folder = os.path.dirname(
                os.path.abspath(pdb)
                )
        os.chdir(folder)
        target = workspace.target_path_clean
        # os.system('ls ???_????_*_*.pdb > input_list')
        # inputs = []
        # with open('input_list', 'r') as f:
            # for line in f:
                # inputs.append(line.strip())
        # for pdb in inputs:

        # Load pose and score function
        pose = pose_from_file(pdb)
        ref = create_score_function('ref2015')
        ref_cst = create_score_function('ref2015')
        ref_cst.set_weight(ScoreType.coordinate_constraint, 1.0)

        # Select chain B for selection
        selector = residue_selector.ChainSelector('B')
        interface_selector_str = \
        '''
        <InterfaceByVector name="interface_selector">
           <Chain chains='A'/>
           <Chain chains='B'/>
        </InterfaceByVector>
        '''
        interface_selector = XmlObjects.static_get_residue_selector(interface_selector_str)
        align_score, patchlength, identity = get_alignment_info(alignment_df, pdb)
        if not designed:
            tf = TaskFactory()

            # Set coordinate constraints for backbone
            coord_cst = constraint_generator.CoordinateConstraintGenerator()
            coord_cst.set_sidechain(False)
            constraints = coord_cst.apply(pose)
            for cst in constraints:
                pose.add_constraint(cst)

            # Initialize fastdesign object
            fastdes = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()

            if args['--buns-penalty']:
                # If buried unsat penalty, use this sfxn
                buns_sfxn = '''
                <ScoreFunction name="sfxn" weights="ref2015" >
                    <Reweight scoretype="approximate_buried_unsat_penalty" weight="5.0" />
                    <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.5" />
                    <Set approximate_buried_unsat_penalty_burial_atomic_depth="4.0" />
                    # Set this to false if you don't know where you might want prolines
                    <Set approximate_buried_unsat_penalty_assume_const_backbone="true" />
                </ScoreFunction>
                '''
                prune_str = '''
                 <PruneBuriedUnsats name="prune"
                 allow_even_trades="false" atomic_depth_probe_radius="2.3" atomic_depth_resolution="0.5" atomic_depth_cutoff="4.5" minimum_hbond_energy="-0.2" />

                '''
                tf.push_back(XmlObjects.static_get_task_operation(prune_str))
                # init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro'\
                        # ' -holes:dalphaball {} -corrections::beta_nov16'.format(dalphaball))
                sfxn = XmlObjects.static_get_score_function(buns_sfxn)
                sfxn.set_weight(ScoreType.coordinate_constraint, 1.0)
                # Set the sfxn
                fastdes.set_scorefxn(sfxn)
            else:
                fastdes.set_scorefxn(ref_cst)

            movemap = MoveMap()
            movemap.set_bb_true_range(pose.chain_begin(2),
                    pose.chain_end(2))
            movemap.set_chi_true_range(pose.chain_begin(2),
                    pose.chain_end(2))
            fastdes.set_movemap(movemap)

            or_selector = residue_selector.OrResidueSelector(selector,
                    interface_selector)
            not_selector = residue_selector.NotResidueSelector(or_selector)
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
                    strlist_to_vector1_str(['GLY']),
                    selector)
            notdesign = operation.OperateOnResidueSubset(no_design,
                    interface_selector)

            if identity > align_threshold:
                favornative = \
                '''
                <FavorNativeResidue name="favornative" bonus="1.5"/>

                '''
                favornative_mvr = XmlObjects.static_get_mover(favornative)
                favornative_mvr.apply(pose)

            if args['--keep-good-rotamers']:
                nopack_selector = utils.list_to_res_selector(nopack)
                good_rots = operation.OperateOnResidueSubset(no_packing,
                        nopack_selector)

            tf.push_back(static)
            tf.push_back(notaa)
            tf.push_back(notdesign)
            tf.push_back(upweight_taskop)
            packertask = tf.create_task_and_apply_taskoperations(pose)

            fastdes.set_task_factory(tf)
            fastdes.apply(pose)
            score = ref(pose)
            # Temp change of PDB filename
            pose.dump_pdb(pdb)

        # Get metrics

        # Align & save (just in case - should not be necessary)
        chainB = pose.split_by_chain(2)
        # pymol_mobile = pymol.cmd.load(pdb, 'mobile')
        # pymol.cmd.align('mobile and chain A', 'target')
        # pymol.cmd.save(pdb, 'mobile')
        # pymol.cmd.delete('mobile')
        flexpep_file = pdb
        flexpep_pose = pose

        # Determine helical propensity
        ss_str = Dssp(chainB).get_dssp_secstruct()
        # secstruct = contiguous_secstruct(ss_str)
        percent_helical = ss_str.count('H') / len(ss_str)

        # Define filters
        buns_all = '''
        <BuriedUnsatHbonds 
            name="Buried Unsat [[-]]"
            report_all_heavy_atom_unsats="true" scorefxn="ref2015"
            cutoff="4" residue_surface_cutoff="20.0"
            ignore_surface_res="true" print_out_info_to_pdb="true"
            dalphaball_sasa="1" probe_radius="1.1" confidence="0"
            only_interface="true" />

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
        ia_mover.apply(flexpep_pose)

        # For delta NPSA, get the two chains
        poseA = utils.pose_get_chain(flexpep_pose, 'A')
        poseB = utils.pose_get_chain(flexpep_pose, 'B')

        # Make filter objects
        buns_all_obj = XmlObjects.static_get_filter(buns_all)
        buns_sc_obj = XmlObjects.static_get_filter(buns_sc)
        npsa_obj = XmlObjects.static_get_filter(npsa)
        exposed_obj = XmlObjects.static_get_filter(exposed_hydrophobics)
        packstat_obj = XmlObjects.static_get_filter(packstat)
        sc_obj = XmlObjects.static_get_filter(sc)
        contact_obj = XmlObjects.create_from_string(contact).get_filter('contact')

        # Calculate delta NPSA
        npsa_complex = npsa_obj.report_sm(flexpep_pose)
        npsa_A = npsa_obj.report_sm(poseA)
        npsa_B = npsa_obj.report_sm(poseB)
        delta_npsa = npsa_complex - npsa_A - npsa_B
        # Calculate NPSA for helix
        npsa_obj.set_residue_selector(selector)
        npsa_score = npsa_obj.report_sm(flexpep_pose)

        buns_all_score = buns_all_obj.report_sm(flexpep_pose)
        buns_sc_score = buns_sc_obj.report_sm(flexpep_pose)
        exposed_score = exposed_obj.report_sm(flexpep_pose)
        packstat_score = packstat_obj.report_sm(flexpep_pose)
        sc_score = sc_obj.report_sm(flexpep_pose)
        contact_score = contact_obj.report_sm(flexpep_pose)

        score = ref(flexpep_pose)
        interface_scorer = interface.InterfaceScore(flexpep_pose)
        interface_score = interface_scorer.apply()
        n_hbonds = interface_scorer.n_hbonds

        # Can't pickle a C++ set, so put it in a Python list
        int_set = []
        for interface_resi in ia_mover.get_interface_set():
            int_set.append(interface_resi)

        row = {'patchman_file': pdb_save,
                'name': os.path.basename(flexpep_file),
                'size': flexpep_pose.size(),
                'pose_score': score,
                'interface_score': interface_score,
                'n_hbonds': n_hbonds,
                'shape_complementarity': sc_score,
                'contact_molecular_surface': contact_score,
                'buns_all': buns_all_score,
                'buns_sc': buns_sc_score,
                # Buried NPSA applies only to helix
                'buried_npsa_helix': npsa_score,
                'delta_buried_npsa': delta_npsa,
                # Exposed hydro. applies to whole protein
                'exposed_hydrophobics': exposed_score,
                'packstat': packstat_score,
                'percent_helical': percent_helical,
                'n_interface_residues': ia_mover.get_num_interface_residues(),
                'complex_sasa': ia_mover.get_complexed_sasa(),
                'delta_sasa': ia_mover.get_interface_delta_sasa(),
                'crossterm_energy': ia_mover.get_crossterm_interface_energy(),
                'interface_packstat': ia_mover.get_interface_packstat(),
                'delta_unsat': ia_mover.get_interface_delta_hbond_unsat(),
                'interface_dG': ia_mover.get_interface_dG(),
                'interfac_residues': int_set,
                'alignment_score': align_score,
                'sequence_identity': identity,
                'patch_length': patchlength,
                'residues_witheld': nopack,
                }
        rowlist.append(row)

    df = pd.DataFrame(rowlist)
    print(df)
    pickle_outdir = os.path.join(workspace.focus_dir, 'scores')
    print('Saving in folder {}'.format(pickle_outdir))
    if not os.path.exists(pickle_outdir):
        os.makedirs(pickle_outdir, exist_ok=True)
    df.to_pickle(os.path.join(pickle_outdir,
        'task_{task}.pkl'.format(task=task_id)))

    # if args['--delete']:
        # os.remove(pdb)

if __name__=='__main__':
    main()
