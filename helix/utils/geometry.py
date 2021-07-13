import numpy as np


def normalize(v):
    '''Normalize a vector based on its 2 norm.'''
    if 0 == np.linalg.norm(v):
        return v
    return v / np.linalg.norm(v)

def angle(v1, v2):
    '''Return the angle between two vectors in radian.'''
    cos = max(-1, min(1, np.dot(normalize(v1), normalize(v2))))
    return np.arccos(cos)

def dihedral(p1, p2, p3, p4):
    '''Return the dihedral defined by 4 points in 
    range [-pi, pi].
    '''
    v1 = normalize(p2 - p1)
    v2 = normalize(p3 - p2)
    v3 = normalize(p4 - p3)

    n1 = normalize(np.cross(v1, v2))
    n2 = normalize(np.cross(v2, v3))

    c = np.dot(n1, n2)
    s = np.dot(v2, np.cross(n1, n2))

    return np.arctan2(s, c) 

def perpendicular_vector(v):
    '''Get a perpendicular vector to a vector v.'''
    v = normalize(v)

    if np.absolute(v[0]) > 0.1:
        return normalize(np.array([v[1], -v[0], 0]))
    else:
        return normalize(np.array([0, v[2], -v[1]]))

def create_frame_from_three_points(p1, p2, p3):
    '''Create a left-handed coordinate frame from 3 points. 
    The p2 is the origin; the y-axis is the vector from p2 to p3; 
    the z-axis is the cross product of the vector from p2 to p1
    and the y-axis.
    
    Return a matrix where the axis vectors are the rows.
    '''
    
    y = normalize(p3 - p2)
    z = normalize(np.cross(p1 - p2, y))
    x = np.cross(y, z)
    return np.array([x, y, z])

def rotation_matrix_to_axis_and_angle(M):
    '''Calculate the axis and angle of a rotation matrix.'''
    u = np.array([M[2][1] - M[1][2],
                  M[0][2] - M[2][0],
                  M[1][0] - M[0][1]])

    sin_theta = np.linalg.norm(u) / 2
    cos_theta = (np.trace(M) - 1) / 2

    return normalize(u), np.arctan2(sin_theta, cos_theta)

def rotation_matrix_from_axis_and_angle(u, theta):
    '''Calculate a rotation matrix from an axis and an angle.'''

    u = normalize(u)
    x = u[0]
    y = u[1]
    z = u[2]
    s = np.sin(theta)
    c = np.cos(theta)

    return np.array([[c + x**2 * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                     [y * x * (1 - c) + z * s, c + y**2 * (1 - c), y * z * (1 - c) - x * s ],
                     [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z**2 * (1 - c) ]])

def cartesian_coord_from_internal_coord(p1, p2, p3, d, theta, tau):
    '''Calculate the cartesian coordinates of an atom from 
    three internal coordinates and three reference points.
    '''
    axis1 = np.cross(p1 - p2, p3 - p2)
    axis2 = p3 - p2

    M1 = rotation_matrix_from_axis_and_angle(axis1, theta - np.pi)
    M2 = rotation_matrix_from_axis_and_angle(axis2, tau)

    return p3 + d * np.dot(M2, np.dot(M1, normalize(p3 - p2)))

def change_angle(p1, p2, p3, theta):
    '''Change the angle between points p1, p2 and p3 to theta.
    Keep p1, p2 and the plane formed by p1, p2 and p3 fixed.
    Return the new position of p3.
    '''
    x = normalize(p1 - p2)
    v = p3 - p2
    y = normalize(v - np.dot(v, x) * x)

    return p2 + np.linalg.norm(v) * (np.cos(theta) * x + np.sin(theta) * y) 

def rotation_matrix_to_superimpose_two_vectors(v1, v2, theta=0):
    '''Get a rotation matrix that superimpose v1 to v2.
    Because there are infinite number of matrices that can do
    so, change the value of theta to get different results.
    '''
    v1 = normalize(v1)
    v2 = normalize(v2)

    axis = np.cross(v1, v2)
    sin_ang = np.linalg.norm(axis)
    cos_ang = np.dot(v1, v2)

    if np.linalg.norm(axis) < 0.01:
        M = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    else:
        M = rotation_matrix_from_axis_and_angle(axis, np.arctan2(sin_ang, cos_ang))

    return np.dot(rotation_matrix_from_axis_and_angle(v2, theta), M)

def point_lists_to_screw_transformation(P1, P2):
    '''Calculate the screw transformation that transforms
    point list P1 to point list P2. Both P1 and P2 have
    3 points respectively.
    '''
    frame1 = create_frame_from_three_points(P1[0], P1[1], P1[2])
    frame2 = create_frame_from_three_points(P2[0], P2[1], P2[2])

    M = np.dot(np.transpose(frame2), frame1)
    t = P2[0] - np.dot(M, P1[0])

    return M, t


def get_screw_transformation(axis, theta, pitch, u):
    '''Get the ration matrix and translate vector of a screw
    transformation from the direction of the screw axis, the 
    screw angle theta, the pitch and a point u on the screw axis.
    If theta > 0, it is a right handed screw and if theta < 0, it is 
    a left handed screw.
    '''
    axis = normalize(axis)
    R = rotation_matrix_from_axis_and_angle(axis, theta)
    t = np.absolute(theta) / (2 * np.pi) * pitch * axis + u - np.dot(R, u)

    return R, t

def get_screw_parameters(M, t):
    '''Get the screw parameters given its Euclidean transformation.
    Return the axis of the screw, the screw angle, the pitch and
    a point on the screw axis. Use the convention that the pitch is
    always positive.
    '''
    axis, theta = rotation_matrix_to_axis_and_angle(M)

    if np.dot(axis, t) < 0:
        axis = -axis
        theta = -theta

    pitch = np.dot(axis, t) * 2 * np.pi / np.absolute(theta)

    u = np.linalg.lstsq(M - np.identity(3), np.dot(axis, t) * axis - t)[0]

    # Let u be perpendicular to the axis

    u = u - np.dot(u, axis) * axis

    return axis, theta, pitch, u

def pitch_angle_to_pitch(pitch_angle, R):
    '''Get the pitch of a screw from its pitch_angle and radius.'''
    return  2 * np.pi * R / np.absolute(np.tan(pitch_angle))

def pitch_to_pitch_angle(pitch, R):
    '''Get the pitch angle of a screw from its pitch and radius.'''
    return np.arctan2(2 * np.pi * R, pitch) 

def get_superimpose_transformation(P1, P2):
    '''Get the superimpose transformation that transfoms a list of
    points P1 to another list of points P2.'''
    if len(P1) != len(P2):
        raise Exception("Sets to be superimposed must have same number of points.")

    com1 = np.mean(P1, axis=0)
    com2 = np.mean(P2, axis=0)

    R = np.dot(np.transpose(np.array(P1) - com1), np.array(P2) - com2)
    V, S, W = np.linalg.svd(R)

    if (np.linalg.det(V) * np.linalg.det(W)) < 0.0:
        V[:, -1] = -V[:, -1]

    M = np.transpose(np.array(np.dot(V, W)))

    return M, com2 - np.dot(M, com1)

def point_line_distance(p, l_p, l_v):
    '''Calculate the distance between a point and a line defined
    by a point and a direction vector.
    '''
    l_v = normalize(l_v)
    u = p - l_p
    return np.linalg.norm(u - np.dot(u, l_v) * l_v)

def point_segment_distance(p, e1, e2):
    '''Calculate the distance between a point and a segment defined
    by two endpoints.
    '''
    if 0 < np.dot(p - e1, e2 - e1) < np.dot(e2 - e1, e2 - e1):
        return point_line_distance(p, e1, e2 - e1)

    return min(np.linalg.norm(p - e1), np.linalg.norm(p - e2))

def intersections_of_circles_on_unit_sphere(v1, v2, theta1, theta2):
    '''Calculate intersection points of two circles on the unit
    sphere. v1, v2 are unit vectors indicating the center of circles.
    theta1 is the anlge between vectors on circle1 and v1, similar 
    to theta2.
    '''
    v1 = normalize(v1)
    v2 = normalize(v2)
    v3 = normalize(np.cross(v1, v2))

    alpha = angle(v1, v2)

    if 0 == np.sin(alpha):
        return None

    x1 = (np.cos(theta1) - np.cos(theta2) * np.cos(alpha)) / np.sin(alpha) ** 2
    x2 = (np.cos(theta2) - np.cos(theta1) * np.cos(alpha)) / np.sin(alpha) ** 2

    p = x1 * v1 + x2 * v2

    if np.linalg.norm(p) > 1:
        return None

    d = np.sqrt(1 - np.dot(p, p)) * v3

    return p + d, p - d
