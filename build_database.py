import os
import pandas as pd
from pyrosetta import init
from pyrosetta import pose_from_file
from scan_helices import PoseScanner

pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
idx = int(os.environ['SGE_TASK_ID']) - 1
print('IDX = {}'.format(idx))
num = 5

def main():
    init('-ignore_unrecognized_res')
    df = pd.DataFrame()
    print('START: {}'.format(idx * num))
    print('STOP: {}'.format(idx * num + num - 1))
    for subdir in sorted(os.listdir(pdb_prefix))[idx*num:idx*num + num - 1]:
        for f in os.listdir(
                os.path.join(
                    pdb_prefix, subdir
                    )
                ):
            if f.endswith('.ent.gz'):
                print('Scanning {}'.format(f))
                path = os.path.join(pdb_prefix, subdir, f)
                pdb = f[3:7]
                try:
                    pose = pose_from_file(path)
                    scanner = PoseScanner(pose)
                    helices = pd.DataFrame(
                            scanner.scan_pose_helices(name=pdb)
                            )

                    pd.concat([df, helices], ignore_index=True)
                except:
                    print("Error scanning {}".format(f))

    df.to_pickle('dataframes/{}.pkl'.format(
        idx
        ))
    df.to_csv('dataframes/{}.csv'.format(
        idx
        ))

if __name__=='__main__':
    main()
