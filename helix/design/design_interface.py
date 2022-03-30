'''
Design the interface and surroundings of a matched scaffold.
Requires score dataframe and match dataframe.
'''
import sys, os
import pandas as pd
from helix import workspace as ws
from pyrosetta import init
from pyrosetta import pose_from_file


def resmap(pose1, pose2, pose1_start, pose1_stop, pose2_start, pose2_stop):
    '''
    Make a map of CA distances for all residues in pose1 to all residues in pose2
    '''
    resmap = {}
    for res1 in range(pose1_start, pose1_stop + 1):
        if res1 not in resmap:
            resmap[res1] = {}
        for res2 in range(pose2_start, pose2_stop + 1):
            if res2 in resmap and res1 in resmap[res2]:
                xyz1 = pose1.residue(res1).xyz('CA')
                xyz2 = pose2.residue(res2).xyz('CA')
                assert(
                    euclidean_distance(xyz1, xyz2) == resmap[res2][res1]
                )
                resmap[res1][res2] = resmap[res2][res1]
            else:
                xyz1 = pose1.residue(res1).xyz('CA')
                xyz2 = pose2.residue(res2).xyz('CA')
                resmap[res1][res2] = euclidean_distance(xyz1, xyz2)

    if len(resmap) > 0:
        resmap = pd.DataFrame(resmap).fillna(0).unstack().reset_index()
        resmap.columns = ['res1', 'res2', 'dist']
    else:
        resmap = None

    return resmap


def prep_design(workspace, export_group):
    '''
    1. Open dataframe, exported scaffold pose, and helix pose
    2. Transfer residues/ROTAMERS over to scaffold pose
        a. For each helix residue, if <0.5 A from scaffold residue, transfer rotamer over with residue type bonus

    Expects a dataframe group where all rows have the same "superimposed_file"
    '''

    for idx, row in export_group.iterrows():
        # print(row['design_file'])
        designpose = pose_from_file(
            os.path.join(
                workspace.root_dir, row.design_file
            )
        )
        print(designpose.chain_sequence(2))


def test_prep():
    df = pd.read_pickle(sys.argv[1])
    init()
    workspace = ws.workspace_from_dir(sys.argv[1])
    for name, group in df.groupby('superimposed_file'):
        print(name)
        prep_design(workspace, group)

if __name__=='__main__':
    test_prep()