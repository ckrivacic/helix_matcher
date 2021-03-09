import geometry
from alpha_helix import helix_direction
from pyrosetta import rosetta
from pyrosetta import get_fa_scorefxn
from pyrosetta import init
from pyrosetta import pose_from_file
import sys
from roseasy.utils import numeric
import numpy as np
import pandas as pd
from statistics import median


def pythonize_vector(v):
    out = []
    for a in v:
        out.append(a)
    return out


def contiguous_secstruct(ss_str):
    ss_positions = {}
    ss_positions['H'] = []
    ss_positions['L'] = []
    ss_positions['E'] = []

    start = 0
    for i in range(0, len(ss_str)):
        if ss_str[i] != ss_str[start]:
            ss_positions[ss_str[start]].append((start + 1, i))
            start = i

        if i + 1== len(ss_str):
            ss_positions[ss_str[start]].append((start + 1, i + 1))

    return ss_positions


def mean(l):
    return sum(l) / len(l)


def find_resis_centroid(resis):
    x = []
    y = []
    z = []
    for resi in resis:
        for atom in resi:
            x.append(resi[atom][0])
            y.append(resi[atom][1])
            z.append(resi[atom][2])

    return np.array([mean(x), mean(y), mean(z)])

def match_2_helices():
    '''
    Just notes for now. To do:
    1. Align one helix with a helix on the scaffold.
    2. From top-down view, find angle between second query helix and
    second scaffold helix.
    3. Do some math to rotate around first query helix. (This is the
    hard part I think)

    Alternatively:
    1. Align 2 query helix vectors with 2 scaffold helix vectors.
    2. Calculate RMSD or something.
    3. From here look for any other matches from same rotation.
    4. Do this with every pair of helices.
    '''
    return


def final_vector(direction, length, centroid):
    # rotation_normalizer = np.cross(
            # np.array([0,0,1]),
            # direction
            # )
    # if np.array_equal(rotation_normalizer, [0,0,0]):
        # rotation_normalizer = np.array([1,0,0])

    # assert np.dot(rotation_normalizer, np.array([0,0,1])) == 0
    # assert np.dot(rotation_normalizer, direction) == 0
            
    # rotation_norm_1 = rotation_normalizer + centroid
    # rotation_norm_2 = centroid - rotation_normalizer

    vector = direction * length
    line_center = vector / 2
    top = vector + centroid - line_center
    bot = centroid - line_center
    # return np.array([bot, rotation_norm_1, rotation_norm_2, top])
    return np.array([bot, top])


def plot_resis(resis, vector):
    """
    Test function that plots a dictionary of residues and a vector that
    describes their orientation
    """
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colordict = {
            'n': 'blue',
            'c': 'black',
            'ca': 'brown'
            }
    for atom in ['n', 'ca', 'c']:
        for resi in resis:
            ax.scatter(resi[atom][0], resi[atom][1], resi[atom][2],
                    color=colordict[atom])

    # centroid = find_resis_centroid(resis)
    # vector_moved = vector + centroid
    x = [point[0] for point in vector]
    y = [point[1] for point in vector]
    z = [point[2] for point in vector]
    # x = [vector[0][0], vector[1][0]]
    # y = [vector[0][1], vector[1][1]]
    # z = [vector[0][2], vector[1][2]]
    ax.plot(x, y, z, color='darkgray', linewidth=4)
    # x = [0, vector[0]]
    # y = [0, vector[1]]
    # z = [0, vector[2]]
    # ax.scatter(x, y, z, color='black')

    plt.show()


class Surface(object):
    def __init__(self):
        # self.sfxn = get_fa_scorefxn()
        # self.neighbor_selector = rosetta.core.select.residue_selector.NumNeighborsSelector(
                    # 2, 5.0
                # )
        self.selector = rosetta.core.select.residue_selector.LayerSelector()
        self.selector.set_layers(False, False, True)


    def apply(self, pose):
        return self.selector.apply(pose)


def find_avg_helix_direction(resis):
    direction_vectors = []
    for i in range(0, len(resis) - 2):
        direction_vectors.append(helix_direction(
            resis[i], resis[i+1], resis[i+2]
            ))
    return np.array(direction_vectors).mean(axis=0)


def helix_length(pose, helix_tuple):
    res1 = pose.residue(helix_tuple[0])
    res2 = pose.residue(helix_tuple[1])

    return numeric.euclidean_distance(res1.xyz('CA'), res2.xyz('CA'))


class PoseScanner(object):
    def __init__(self, pose):
        self.selector = Surface()
        self.pose = pose


    def scan_pose_helices(self, name=None, test=False,
            split_chains=True, path=None):
        """
        Scan a pose for helices, then find their centroid, direction,
        length, name, and solvent accessibility.
        To do: solvent accessbility, name (can maybe be done outside this
        function).
        """
        if not name:
            name = self.pose.pdb_info().name()
            print('NAME')
            print(name)
        if split_chains:
            chains = self.pose.split_by_chain()
        else:
            chains = [self.pose]
        # Calculate which residues in pose are at surface only once
        ch = 1
        helices_found = []
        for pose in chains:
            dssp = rosetta.core.scoring.dssp.Dssp(pose)
            positions = contiguous_secstruct(dssp.get_dssp_secstruct())
            surface = np.array(
                    pythonize_vector(self.selector.apply(pose))
                    )
            for helix in positions['H']:
                # print(helix)
                if helix[1] - helix[0] > 1:
                    helix_info = {}
                    helix_info['start'] = helix[0]
                    helix_info['stop'] = helix[1]
                    resis = []
                    positions = np.arange(helix[0], helix[1] + 1)
                    # med = int(median(positions)) - 1
                    # start = min(helix[1] - 3, med)
                    for i in range(helix[0], helix[1] + 1):
                        res = {}
                        for atom in ['n', 'ca', 'c']:
                            res[atom] = numeric.xyzV_to_np_array(pose.residue(i).xyz(atom.upper()))
                        resis.append(res)

                    avg_direction = find_avg_helix_direction(resis)
                    # This gives us direction only; need to also get center of
                    # mass and maybe length of helix to fully characterize position
                    # helix_info['direction'] = helix_direction(resis[0], resis[1], resis[2])
                    helix_info['direction'] = avg_direction
                    helix_info['centroid'] = find_resis_centroid(resis)
                    helix_info['nres'] = helix[1] - helix[0]
                    helix_info['length'] = helix_length(pose, helix)
                    if split_chains:
                        helix_info['name'] = name + '_' + str(ch)
                    else:
                        helix_info['name'] = name
                    helix_info['vector'] = final_vector(helix_info['direction'], 
                            helix_info['length'], helix_info['centroid'])
                    helix_info['surface'] = surface[helix[0]-1:helix[1]-1]
                    helix_info['percent_exposed'] =\
                            np.count_nonzero(helix_info['surface']) /\
                            len(helix_info['surface'])
                    helix_info['chain'] = pose.pdb_info().pose2pdb(helix_info['start'])
                    if path:
                        helix_info['path'] = path
                    helices_found.append(helix_info)
                    if test:
                        plot_resis(resis, helix_info['vector'])
            ch += 1

        return helices_found


def main():
    init()
    pose = pose_from_file('test_files/6r9d.cif')
    # pose = pose_from_file('test_cas9.pdb.gz')
    import time
    start = time.time()
    scanner = PoseScanner(pose)
    helices = scanner.scan_pose_helices(test=False)
    print('CALC TIME = {}'.format(
        time.time() - start
        ))
    helices = pd.DataFrame(helices)
    print(helices)
    print(helices['name'])
    helices.to_pickle('out.pkl')
    helices.to_csv('out.csv')

if __name__=='__main__':
    main()
