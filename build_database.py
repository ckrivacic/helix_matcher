import os
from pyrosetta import init
from pyrosetta import pose_from_file
from scan_helices import PoseScanner

pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'

def main():
    df = pd.DataFrame()
    for root, dirs, files in os.walk(pdb_prefix):
        for f in files:
            if f.endswith('.ent.gz'):
                path = os.path.join(root, f)
                pdb = f[3:7]
                pose = pose_from_file(path)
                scanner = PoseScanner(pose)
                helices = pd.DataFrame(
                        scanner.scan_pose_helices(name=pdb)
                        )

                pd.concat([df, helices], ignore_index=True)
