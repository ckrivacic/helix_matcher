import os
from pyrosetta import init
from pyrosetta import pose_from_file
from scan_helices import PoseScanner

pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
idx = 0

def main():
    init('-ignore_unrecognized_res')
    df = pd.DataFrame()
    for subdir in os.listdir(pdb_prefix)[idx*100:idx*100 + 100 - 1]:
        for f in os.listdir(
                os.path.join(
                    pdb_prefix, subdir
                    )
                ):
            if f.endswith('.ent.gz'):
                path = os.path.join(root, f)
                pdb = f[3:7]
                pose = pose_from_file(path)
                scanner = PoseScanner(pose)
                helices = pd.DataFrame(
                        scanner.scan_pose_helices(name=pdb)
                        )

                pd.concat([df, helices], ignore_index=True)

    pd.to_pickle('dataframes/{}.pkl'.format(
        idx
        ))
    pd.to_csv('dataframes/{}.csv'.format(
        idx
        ))
