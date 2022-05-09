"""
Script to analyze all structures from the Baker lab design method.

Usage:
    analyze_structures.py <folder> [options]

Options:
    --task=INT  Only run a certain task
    --designs-per-task=INT  How many designs to analyze per task  [default: 1]
    --benchmark  Run analyhsis on benchamrk pdbs
"""

from pyrosetta import rosetta
# from pyrosetta.rosetta.core.pack.task import TaskFactory
# from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.select import residue_selector
from pyrosetta.rosetta.core.scoring.dssp import Dssp
# from pyrosetta.rosetta.core.kinematics import MoveMap
# from pyrosetta.rosetta.core.scoring import ScoreType
# from pyrosetta.rosetta.protocols import constraint_generator
from pyrosetta import init
from pyrosetta import pose_from_file
from pyrosetta import create_score_function
import pyrosetta
from helix.patchman.pythondesign import SpecialRotDesign
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


def main():
    args = docopt.docopt(__doc__)

    if args['--task']:
        task_id = int(args['--task']) - 1
    elif 'task' in os.environ:
        task_id = int(os.environ['task']) - 1
    else:
        task_id = 0
    rosetta_dir = os.path.expanduser('~/software/rosetta/')
    dalphaball = os.path.join(rosetta_dir,
                              'source', 'external', 'DAlpahBall',
                              'DAlphaBall.gcc')
    init('-total_threads 1 -ex1 -ex2 -use_input_sc -ex1aro' \
         ' -holes:dalphaball {} -ignore_unrecognized_res -detect_disulf false'.format(dalphaball))

    inputs = sorted(glob.glob(args['<folder>'] + '/*.pdb.gz'))
    nstruct = int(args['--designs-per-task'])
    total_jobs = len(inputs)
    print('TOTAL JOBS: {}'.format(total_jobs), flush=True)

    start = task_id * nstruct
    stop = task_id * nstruct + nstruct
    if stop > len(inputs) - 1:
        stop = len(inputs) - 1
    if stop < start:
        sys.exit('Nothing to do here')
    print(start, flush=True)
    print(inputs[start], flush=True)

    rowlist = []
    for input_idx in range(start, stop):
        pdb = inputs[input_idx]
        pose = pose_from_file(pdb)
        if pose.num_chains() < 2:
            continue
        row = analyze_pose(pose, chA='A', chB='B')
        rowlist.append(row)

    df = pd.DataFrame(rowlist)
    pickle_outdir = args['<folder>']
    print(df, flush=True)
    print('Saving in folder {}'.format(pickle_outdir), flush=True)
    if not os.path.exists(pickle_outdir):
        os.makedirs(pickle_outdir, exist_ok=True)
    dataframe_out = os.path.join(pickle_outdir, f'task_{task_id}.pkl')
    df.to_pickle(dataframe_out)

def analyze_pose(pose, chA='A', chB='B'):
    flexpep_pose = pose.clone()
    ref = create_score_function('ref2015')
    ref(flexpep_pose)

    # Need to clear sequence constraints to split the pose
    pose.remove_constraints()
    pose.clear_sequence_constraints()
    # chainB = pose.split_by_chain(2)
    chainB = utils.pose_get_chain(pose, chB)
    chainA = utils.pose_get_chain(pose, chA)
    pose = utils.pose_get_chain(pose, chA)
    ref(chainA)
    ref(chainB)
    # This way we only have chain A and chain B in cases where we are looking at a multimer (particularly in the
    # case of benchamrk examples)
    rosetta.core.pose.append_pose_to_pose(chainA, chainB)
    print('COMBINED POSE CHAINS')
    print(pose.num_chains())
    ref(pose)

    # Determine helical propensity
    ss_str = Dssp(chainB).get_dssp_secstruct()
    # secstruct = contiguous_secstruct(ss_str)
    percent_helical = ss_str.count('H') / len(ss_str)
    print("SETTING UP FILTERS", flush=True)

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
        f'''
    <RESIDUE_SELECTORS>
        <Chain name="chA" chains="{chA}"/>
        <Chain name="chB" chains="{chB}"/>
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
    ia_mover.apply(flexpep_pose)

    # For delta NPSA, get the two chains
    poseA = chainA
    poseB = chainB
    ref(poseA)
    ref(poseB)

    # Make filter objects
    buns_all_obj = XmlObjects.static_get_filter(buns_all)
    buns_sc_obj = XmlObjects.static_get_filter(buns_sc)
    npsa_obj = XmlObjects.static_get_filter(npsa)
    exposed_obj = XmlObjects.static_get_filter(exposed_hydrophobics)
    packstat_obj = XmlObjects.static_get_filter(packstat)
    sc_obj = XmlObjects.static_get_filter(sc)
    contact_obj = XmlObjects.create_from_string(contact).get_filter('contact')

    selector = residue_selector.ChainSelector(chB)

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

    row = {
           'file': pdb,
            'folder': os.path.dirname(pdb),
           'name': os.path.basename(pdb),
           # 'protocol': designtype,
            'protocol': os.path.basename(args['<folder>']),
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
           }
    return row


if __name__=='__main__':
    main()