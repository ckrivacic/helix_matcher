import os, sys
import pandas as pd
from pyrosetta import init
from pyrosetta import pose_from_file
from scan_helices import PoseScanner
import gzip

pdb_prefix = '/wynton/home/database/pdb/remediated/pdb/'
idx = int(os.environ['SGE_TASK_ID']) - 1
print('IDX = {}'.format(idx))
num = 2000


def test_iterate():
    f = gzip.open('test_files/nrpdb.gz', 'rb')
    lines = f.readlines()
    for line in lines:
        line = line.decode('utf-8')
        if not line.startswith('#'):
            get_pose(str(line))

def get_pose(line):
    fields = list(filter(None, line.split(' ')))
    pdb = fields[0].lower()
    chain = fields[1]
    rep = fields[5]
    if rep:
        path = os.path.join(pdb_prefix, pdb[1:3], 'pdb{}.ent.gz'.format(
            pdb
            ))
        pose = pose_from_file(path)
        pdbinfo = pose.pdb_info()
        i = 1
        rosetta_number = pdbinfo.pdb2pose(
            chain, i)
        while rosetta_number == 0:
            i += 1
            rosetta_number = pdbinfo.pdb2pose(
                    chain, i
                    )
            if i > pose.size():
                break
        chain = pose.chain(rosetta_number)

        return pose.split_by_chain(chain)
    else:
        return None

def main():
    init('-ignore_unrecognized_res')
    df = pd.DataFrame()
    start = idx * num
    stop = idx * num + num - 1
    print('START: {}'.format(start))
    print('STOP: {}'.format(stop))
    with gzip.open('test_files/nrpdb.gz', 'rb') as f:
        lines = f.readlines()[start:stop]
    errors = []
    for line in lines:
        line = line.decode('utf-8')
        if not line.startswith('#'):
            try:
                print('Opening from line {}'.format(line))
                sys.stdout.flush()
                pose = get_pose(str(line))
                if pose:
                    scanner = PoseScanner(pose)
                    helices = pd.DataFrame(
                            scanner.scan_pose_helices(name=pdb)
                            )

                    df = pd.concat([df, helices], ignore_index=True)
            except Exception as e:
                print("Error scanning line: \n{}".format(line))
                print('Error was:')
                print(e)
                sys.stdout.flush()
                errors.append(line)
    # for subdir in sorted(os.listdir(pdb_prefix))[idx*num:idx*num + num - 1]:
        # for f in os.listdir(
                # os.path.join(
                    # pdb_prefix, subdir
                    # )
                # ):
            # if f.endswith('.ent.gz'):
                # print('Scanning {}'.format(f))
                # path = os.path.join(pdb_prefix, subdir, f)
                # pdb = f[3:7]
                # try:
                    # pose = pose_from_file(path)
                    # scanner = PoseScanner(pose)
                    # helices = pd.DataFrame(
                            # scanner.scan_pose_helices(name=pdb)
                            # )

                    # df = pd.concat([df, helices], ignore_index=True)
                # except:
                    # print("Error scanning {}".format(f))

    df.to_pickle('dataframes/{}.pkl'.format(
        idx
        ))
    df.to_csv('dataframes/{}.csv'.format(
        idx
        ))

    errorlog = os.path.join('dataframes', 'errors',
            'helix_scan.e{}'.format(idx))
    with open(errorlog, 'w') as f:
        for err in errors:
            f.write(err + '\n')


if __name__=='__main__':
    main()
    # test_iterate()
