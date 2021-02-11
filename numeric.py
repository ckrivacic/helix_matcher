import numpy as np
import math
import pyrosetta
from pyrosetta import rosetta


class Transformation(object):
    """Class for calculating and storing rotation & transformation information"""

    def __init__(self, xyz1, xyz2, centroid=None):
        self.xyz1 = xyz1
        self.xyz2 = xyz2
        self.rotation, self.translation = \
                self.get_superimpose_transformation( 
                self.xyz1, self.xyz2, centroid=centroid)

    def get_superimpose_transformation(self, P1, P2, centroid=None):
        '''Get the superimpose transformation that transfoms a list of
        points P1 to another list of points P2.
        From XingJie Pan
        '''
        if len(P1) != len(P2):
            raise Exception("Sets to be superimposed must have same number of points.")

        com1 = np.mean(P1, axis=0)
        com2 = np.mean(P2, axis=0)

        if centroid is None:
            centroid = com1

        R = np.dot(np.transpose(np.array(P1) - centroid), np.array(P2) - com2)
        V, S, W = np.linalg.svd(R)

        if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
            V[:, -1] = -V[:, -1]

        M = np.transpose(np.array(np.dot(V, W)))

        return M, com2 - np.dot(M, com1)

    def apply(self, vector):
        rotated = np.transpose(np.dot(self.rotation,
            np.transpose(vector)))
        for i in range(0, len(rotated)):
            rotated[i] = rotated[i] + self.translation
        return rotated


def centroid(vector1, vector2):
    '''
    Finds centroid of two sets of points
    '''
    cen1 = np.mean(vector1, axis=0)
    cen2 = np.mean(vector2, axis=0)
    return np.mean([cen1, cen2], axis=0)


def angle(a, b, c):
    '''
    Finds the angle between 3 points in degrees.
    '''
    ba = a - b
    bc = c - b

    cos = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cos)
    return np.degrees(angle)



def dihedral(a, b, c, d):
    '''
    Finds the dihedral angle between 4 points in degrees.
    '''
    b0 = -1.0 * (b - a)
    b1 = c - b
    b2 = d - c

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1 = b0 minus
    # component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1 = b2 minus
    # component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))



def line_angle(line1, line2):
    '''
    Finds the angle between two lines, each consisting of two points.
    E.g. line1 = [[1,2,3],[5,6,7]]
    '''
    # Move to origin
    vector1 = line1[1] - line1[0]
    vector2 = line2[1] - line2[0]
    # Normalize
    vector1 = vector1 / np.linalg.norm(vector1)
    vector2 = vector2 / np.linalg.norm(vector2)
    return np.arccos(np.clip(np.dot(vector1, vector2), -1.0, 1.0))


def apply_transformation(Transformation, template_coordinate_set):
    return np.dot(template_coordinate_set, Transformation.rotation.T) +\
            Transformation.translation

def xyzV_to_np_array(xyz):
    return np.array([xyz.x, xyz.y, xyz.z])

def np_array_to_xyzV(a):
    return rosetta.numeric.xyzVector_double_t(a[0], a[1], a[2])

def intlist_to_vector1_size(a):
    vector = rosetta.utility.vector1_unsigned_long()
    for item in a:
        vector.append(item)
    return vector

def xyzM_to_np_array(M):
    return np.array([[M.xx, M.xy, M.xz],
                     [M.yx, M.yy, M.yz],
                     [M.zx, M.zy, M.zz]])

def np_array_to_xyzM(a):
    return rosetta.numeric.xyzMatrix_double_t.rows(
            a[0][0], a[0][1], a[0][2],
            a[1][0], a[1][1], a[1][2],
            a[2][0], a[2][1], a[2][2])

def mult_np_transformation(T1, T2):
    '''Multiply two numpy rigid body transformations'''
    M1, v1 = T1
    M2, v2 = T2
    
    return np.dot(M1, M2), np.dot(M1, v2) + v1

def inverse_np_transformation(T):
    '''Inverse an numpy rigid body transformation.'''
    M, v = T
    
    invM = np.linalg.inv(M)   
    return invM, - np.dot(invM, v)


def RMSD(points1, poinsts2):
    '''Calcualte RMSD between two lists of numpy points.'''
    diff = [points1[i] - poinsts2[i] for i in range(len(points1))]
    return np.sqrt(sum(np.dot(d, d) for d in diff) / len(diff))


def xyz_from_3d_array(array):
    """
    Takes a 3-dimensional numpy array and returns lists of x, y, and z
    coordinates.
    """
    x = array[:,0]
    y = array[:,1]
    z = array[:,2]

    return x,y,z


def xyz_to_array(xyz):
    """
    Convert a list of strings representing a 3D coordinate to floats and return
    the coordinate as a ``numpy`` array.
    """
    return np.array([float(x) for x in xyz])


def euclidean_distance(xyz1, xyz2):
    """
    Simple function for calculating euclidean distance between two points.
    """
    dist = [(a - b)**2 for a,b in zip(xyz1, xyz2)]
    return math.sqrt(sum(dist))


def backbone_rmsd(rotamer, residue,
        alignment_atoms):
    """
    Measure backbone RMSD between a rotamer and the nearest residue on
    the design protein.
    """

    distances = np.array([])
    for atom in alignment_atoms:
        rot_xyz = xyzV_to_np_array(rotamer.xyz(atom))
        near_xyz = xyzV_to_np_array(residue.xyz(atom))
        distances = np.append(distances,euclidean_distance(rot_xyz,near_xyz))

    return np.sqrt((distances**2).mean())
