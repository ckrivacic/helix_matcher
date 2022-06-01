'''
Make a PyMOL session for figure 4A, showing

Usage:
    fig_4a.py <match_workspace> <model> [options]

Options:
'''

import helix.workspace as ws
from helix.utils import utils
from helix.commands.plot_helix_distribution import get_json
import pymol
import os
import sys
import docopt
from pyrosetta import init
from pyrosetta import pose_from_file
from helix.utils.colors import pymol_palette as palette
from scipy.spatial import KDTree


def selstr_from_reslist(reslist, name, chain=None, chains=None):
    selstr = f"{name} and ("
    for idx, resnum in enumerate(reslist):
        selstr += f"(resi {resnum}"
        if chain:
            selstr += f" and chain {chain}) or "
        elif chains:
            selstr += f" and chain {chains[idx]}) or "
        else:
            selstr += ") or "
    selstr = selstr[:-4]
    selstr += ")"
    return selstr


def show_interface_hbonds(model_name_1, model_name_2, chain_1, chain_2):
    '''Show hbonds between two chains.'''
    pymol.cmd.select("prot1", f"{model_name_1} and chain {chain_1}")
    pymol.cmd.select("prot2", f"{model_name_2} and chain {chain_2}")
    pymol.cmd.h_add('prot1')
    pymol.cmd.h_add('prot2')
    pymol.cmd.select('don', '(elem n,o and (neighbor hydro))')
    pymol.cmd.select('acc', '(elem o or (elem n and not (neighbor hydro)))')
    pymol.cmd.distance('HBA', '(prot1 and acc)', '(prot2 and don)', 3.2)
    pymol.cmd.distance('HBD', '(prot1 and don)', '(prot2 and acc)', 3.2)
    pymol.cmd.delete('don')
    pymol.cmd.delete('acc')
    pymol.cmd.hide('everything', 'hydro')
    pymol.cmd.hide('labels', 'HBA')
    pymol.cmd.hide('labels', 'HBD')


def load_patchman_ex(workspace, row):
    design_file = os.path.join(workspace.root_dir, row['design_file'])
    design_name = os.path.basename(row['design_file']).split('.')[0]

    name_split = design_name.split('_')
    parent_pdb = name_split[1]
    stretch_start = int(name_split[3])
    stretch_stop = int(name_split[4])

    pymol.cmd.fetch(parent_pdb)
    parent_pdb_path = os.path.join(os.path.expanduser(
        f'~/fetch/{parent_pdb}.cif'
    ))
    pose = pose_from_file(parent_pdb_path)
    pdbinfo = pose.pdb_info()

    selstr = f"{parent_pdb} and ("
    for resnum in range(stretch_start, stretch_stop + 1):
        res_info = pdbinfo.pose2pdb(resnum).split()
        selstr += f"(resi {res_info[0]} and chain {res_info[1]}) or "
    selstr = selstr[:-4]
    selstr += ")"
    pymol.cmd.color(palette['lightgray'], f"{parent_pdb}")

    return selstr, parent_pdb


def make_session(workspace, model_path, input_path, patchman_rows):
    insertions = get_json(workspace, input_path)

    model_name = os.path.basename(model_path).split('.')[0]
    input_name = os.path.basename(input_path).split('.')[0]
    pymol.cmd.load(model_path, model_name)
    pymol.cmd.color(palette['blue'], f"{model_name} and chain A and name c*")
    pymol.cmd.color(palette['white'], f"{model_name} and chain B and name c*")
    pymol.cmd.load(input_path, input_name)
    pymol.cmd.align(f"{model_name} and chain B", f"{input_name} and chain B")
    pymol.cmd.hide('everything', f"{input_name} and chain B")
    pymol.cmd.color(palette['blue'], f"{input_name} and chain A and name c*")
    pymol.cmd.color(palette['white'], f"{input_name} and chain B and name c*")

    pymol.cmd.create("model_surface", f"({model_name}) and chain A")
    pymol.cmd.create("target_surface", f"({model_name}) and chain B")
    pymol.cmd.color(palette['blue'], 'model_surface')
    pymol.cmd.color(palette['white'], 'target_surface')

    for insertion in insertions:
        start_res = insertion['start']
        stop_res = insertion['stop']
        pymol.cmd.color(palette['teal'], f"({model_name} or {input_name}) and resi {start_res}-{stop_res} and name c*")
        pymol.cmd.color(palette['teal'], f"model_surface and resi {start_res}-{stop_res}")

    pymol.cmd.load(workspace.target_path, 'target')
    pymol.cmd.color(palette['white'], 'target')

    pymol.cmd.center(f"{model_name} and chain A")


    helices = []
    for idx, row in patchman_rows.iterrows():
        pymol_space = {'patch_resnums': [], 'patch_xyz': [], 'pdb_resnums': [], 'pdb_xyz': [], 'pdb_chains': []}
        design_file = os.path.join(workspace.root_dir, row['design_file'])
        design_name = os.path.basename(row['design_file']).split('.')[0]
        helices.append(design_name)
        pymol.cmd.load(design_file, design_name)
        pymol.cmd.align(f"{design_name} and chain A", f"{input_name} and chain B")
        pymol.cmd.hide('everything', f"{design_name} and chain A")
        pymol.cmd.color(palette['orange'], f"{design_name} and name c* and chain B")

        patch_dir = os.path.dirname(os.path.dirname(design_file))
        patch_number = os.path.basename(os.path.dirname(patch_dir)).split('_')[1]
        patch_file = os.path.join(patch_dir, f"{patch_number}_target.pdb")
        patch_name = os.path.basename(patch_file).split('.')[0]
        pymol.cmd.load(patch_file, patch_name)
        pymol.cmd.color(palette['pink'], patch_name)

        pymol.cmd.iterate(f"{patch_name} and name CA", "patch_resnums.append(resi)", space=pymol_space)
        pymol.cmd.iterate_state(1, f"{patch_name} and name CA", "patch_xyz.append([x, y, z])", space=pymol_space)
        pymol.cmd.select(f"{patch_number}_patch", selstr_from_reslist(pymol_space['patch_resnums'], 'target'))
        pymol.cmd.color(palette['pink'], f"{patch_number}_patch and name c*")

        selection, pdbid = load_patchman_ex(workspace, row)
        pymol.cmd.select(f"patch_{patch_number}_stretch", selection)
        pymol.cmd.super(f"patch_{patch_number}_stretch", f"{design_name} and chain B", cycles=0)
        pymol.cmd.color(palette['orange'], f"patch_{patch_number}_stretch and name c*")

        pymol.cmd.iterate_state(1, f"{pdbid} and name CA", "pdb_resnums.append(resi)", space=pymol_space)
        pymol.cmd.iterate_state(1, f"{pdbid} and name CA", "pdb_chains.append(chain)", space=pymol_space)
        pymol.cmd.iterate_state(1, f"{pdbid} and name CA", "pdb_xyz.append([x, y, z])", space=pymol_space, atomic=0)
        tree = KDTree(pymol_space['pdb_xyz'])

        match_resnums = []
        match_chains = []
        for xyz in pymol_space['patch_xyz']:
            nearest = tree.query(xyz)
            index = nearest[1]
            match_resnums.append(pymol_space['pdb_resnums'][index])
            match_chains.append(pymol_space['pdb_chains'][index])
        pymol.cmd.select(f"{patch_number}_match", selstr_from_reslist(match_resnums, pdbid, chains=match_chains))
        pymol.cmd.color(palette['pink'], f"{patch_number}_match and name c*")

        pymol.cmd.set('cartoon_transparency', 0.5, f'{pdbid} and not (*_stretch or *_match)')

    pymol.cmd.center(f"{model_name} and chain A")
    pymol.cmd.hide('everything', 'resn hoh or resn so4 or resn edo or resn cl or resn na or resn mn or resn zn or resn mg or resn ca')
    show_interface_hbonds(model_name, model_name, 'A', 'B')
    pymol.cmd.set('transparency', 0.5)
    pymol.cmd.disable('all')
    pymol.cmd.enable('target')
    for helix in helices:
        pymol.cmd.enable(helix)

    interface_path = '/Users/codykrivacic/software/anaconda3/envs/proteindesign/lib/python3.7/site-packages/pmg_tk/startup/interfaceFinder.py'
    sys.path.insert(0, os.path.dirname(interface_path))
    import interfaceFinder
    interfaceFinder.interfaceResidues(model_name)

    pymol.cmd.show('sticks', "interface and not hydro")
    pymol.cmd.show('surface', 'model_surface')
    pymol.cmd.show('surface', 'target_surface')
    target_name = workspace.basename(workspace._initial_target_path)
    save_name = f'fig4_sessions/{target_name}_{model_name}.pse'
    pymol.cmd.save(save_name)
    os.system(f'pymol {save_name}')


def main():
    if not os.path.exists(os.path.expanduser('~/fetch')):
        os.makedirs(os.path.expanduser('~/fetch'), exist_ok=True)
    pymol.cmd.set('fetch_path', pymol.cmd.exp_path('~/fetch'))
    init()

    args = docopt.docopt(__doc__)
    workspace = ws.workspace_from_dir(args['<match_workspace>'])
    model = args['<model>']
    if len(model.split('/')) > 1:
        model_path = model
    else:
        if not model.startswith('model_'):
            model = 'model_' + model
        if not model.endswith('.pdb.gz'):
            model += '.pdb.gz'
        model_path = os.path.join(workspace.design_dir, model)

    input_model = '_'.join(os.path.basename(model_path).split('_')[:3]) + '.pdb.gz'
    input_path = os.path.abspath(os.path.join(workspace.complex_dir, input_model))
    input_relpath = os.path.relpath(input_path, workspace.root_dir)

    exported_pdbs = utils.safe_load(os.path.join(
        workspace.complex_dir, 'exported_pdbs.pkl'
    ))

    rows = exported_pdbs[exported_pdbs['superimposed_file'] == input_relpath]
    make_session(workspace, model_path, input_path, rows)


if __name__ == '__main__':
    main()