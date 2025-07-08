import numpy as np
import random

from typing import Dict, List, Tuple, Union, Optional
from scipy.spatial import KDTree


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    rad = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    deg = rad*180./np.pi
    if deg > 90.0 :
        deg = 180.0 - deg
    return deg


def quaternion_to_rotation_matrix(Q:List[float]) -> np.ndarray:
    """Covert a quaternion into a full three-dimensional rotation matrix.

    Args:
        Q (List[float]): A 4 element array representing the quaternion (q0,q1,q2,q3)

    Returns:
        np.ndarray: A 3x3 element matrix representing the full 3D rotation matrix. 
            This rotation matrix converts a point in the local reference frame 
            to a point in the global reference frame.
    """

    # Extract the values from Q
    q0 = Q[0]
    q1 = Q[1]
    q2 = Q[2]
    q3 = Q[3]
     
    # First row of the rotation matrix
    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)
     
    # Second row of the rotation matrix
    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)
     
    # Third row of the rotation matrix
    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1
     
    # 3x3 rotation matrix, R
    R = np.array([[r00, r01, r02],
                           [r10, r11, r12],
                           [r20, r21, r22]])
                            
    return R


def random_translation_vector(lower:float=0.0, upper:float=1.0) -> np.ndarray:
    """Generate a random translational vector.

    Args:
        lower (float) : lower bound. Defaults to 0.0
        upper (float) : upper bound. Defaults to 1.0

    Returns:
        np.ndarray: random (3,) vector
    """
    rng = np.random.default_rng()
    return rng.uniform(lower, upper, size=3)



def random_rotation_matrix() -> np.ndarray:
    """Generate a random rotational matrix.

    Returns:
        ndarray : random (3,3) matrix
    """
    axis = np.random.randn(3)
    axis /= np.linalg.norm(axis)

    angle = random.uniform(0, 2 * np.pi)
    # Build the rotation matrix using Rodrigues' formula
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * np.dot(K, K)

    return R


def transform(P:np.ndarray,
              centroid:np.ndarray = np.array([0., 0., 0.]),
              rot:np.ndarray = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]), 
              trans:np.ndarray = np.array([0., 0., 0.])) -> np.ndarray:
    """Rotate and translate input coordinates.

    Args:
        P (np.ndarray): input coordinates.
        centroid (np.ndarray, optional): centroid. Defaults to np.array([0., 0., 0.])
        rot (np.ndarray, optional): rotation matrix. Defaults to np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]).
        trans (np.ndarray, optional): translation vector. Defaults to np.array([0., 0., 0.]).

    Returns:
        np.ndarray: transformed output coordinates.
    """
    Q = np.dot(P-centroid, rot.T) + centroid + trans
    return Q


def kabsch_algorithm(P:np.ndarray, Q:np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Computes the optimal rotation and translation to align two sets of points (P -> Q),
    and their RMSD.

    Examples:
        >>> P = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        >>> Q = np.array([[1, 1, 0], [2, 1, 0], [1, 2, 0]])
        >>> # Q is translated by (1,1,0) 
        >>> (rot, trans, rmsd) = kabsch_transform(P, Q)
        >>> transform(P, rot, trans)

        >>> P = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        >>> Q = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        >>> # Q is a 90-degree rotation of P around the z-axis
        >>> (rot, trans, rmsd) = kabsch_transform(P, Q)
        >>> transform(P, rot, trans)

    Args:
        P (np.ndarray): subject coordinates. Not modified.
        Q (np.ndarray): target coordinates. Not modified.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray, float]: A tuple containing the optimal rotation matrix, the optimal
            translation vector, the centroid of P, and the RMSD.
    """
    
    assert P.shape == Q.shape, "Matrix dimensions must match"

    # Compute centroids
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)

    # Optimal translation
    t = centroid_Q - centroid_P

    # Center the points
    p = P - centroid_P
    q = Q - centroid_Q

    # Compute the covariance matrix
    H = np.dot(p.T, q)

    # SVD
    U, S, Vt = np.linalg.svd(H)
    V = Vt.T
    d = np.linalg.det(np.dot(V, U.T))
    e = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, d]])
    
    # Optimal rotation
    R = np.dot(np.dot(V, e), U.T)

    # RMSD
    rmsd = np.sqrt(np.sum(np.square(np.dot(p, R.T) - q)) / P.shape[0])

    return (R, t, centroid_P, rmsd)


def test_kabsch_algorithm(trials:int=100, rtol:float=1e-5, atol:float=1e-8) -> None:
    results = []
    for i in range(trials):
        rot = random_rotation_matrix()
        trans = random_translation_vector(0.0, 50.0)
        # generate a 10x3 array with random floats between 0 and 1
        P = np.random.rand(10, 3) * 100.0
        Q = transform(P, rot, trans)
        (rot_, trans_, centroid_, rmsd_) = kabsch_algorithm(P, Q)
        res = np.allclose(rot, rot_, rtol=rtol, atol=atol) and np.allclose(trans, trans_, rtol=rtol, atol=atol)
        results.append(res)
    if all(results):
        print("pass")
    else:
        print("failed")



def pca(X):
    # Data matrix X
    n, m = X.shape
    centroid = X.mean(axis=0)
    X = X - centroid # make it 0-centered
    assert np.allclose(X.mean(axis=0), np.zeros(m), atol=1e-5)
    # Compute covariance matrix
    C = np.dot(X.T, X) / (n-1)
    # Eigen decomposition
    eigen_vals, eigen_vecs = np.linalg.eig(C)
    # Project X onto PC space
    # X_pca = np.dot(X, eigen_vecs)
    return eigen_vecs



def hexagonal_grid() -> np.ndarray:
    """Generate XY coordinates of hexagonal grid for membrane MD simulations.

    Returns:
        np.ndarray: (x,y) coordinates
    """
    theta = np.radians(60)
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])

    d = 50.0
    grid = {0: [ np.array([0,0]) ] }
    for shell in [1, 2]:
        v = np.array([0.0, d*shell])
        grid[shell] = [v]
        for i in range(0, 5):
            v = np.dot(R, v)
            if shell > 1:
                delta = (grid[shell][-1] - v)/shell
                for j in range(1, shell):
                    grid[shell].append(j*delta + v)
            grid[shell].append(v)
        if shell > 1:
            v = np.dot(R, v)
            delta = (grid[shell][-1] - v)/shell
            for j in range(1, shell):
                grid[shell].append(j*delta + v)
        grid[shell].append(v)
        
    xy = []
    for shell in grid:
        for item in grid[shell] :
            x_, y_ = item
            xy.append([x_, y_])
    return np.array(xy)



def get_closest_distance(fixed:np.ndarray, 
                         movable:np.ndarray) -> float:
    """Get the closest distance between two sets of coordinates.

    Args:
        fixed (np.ndarray): input coordinates.
        movable (np.ndarray): input coordinates.

    Returns:
        float: the closest distance
    """
    kd = KDTree(fixed)
    distances, indices = kd.query(movable)
    return np.min(distances)



def spherical_points(r:float, rstep:float, azstep:float) -> np.ndarray:
    """3D coordinates evenly distributed on a sphere.

    The spherical coordinates are defined by a radius, polar angle (theta), and azimuthal angle (phi).
    Cartesian coordinates are defined by:
        x = r sin(theta) cos(phi)
        y = r sin(theat) sin(phi)
        z = r cos(theta)
    Because dA = r^2 sin(theta) dtheta dphi = r d(r cos(theta)) dphi = r dz dphi,
    Evenly distributed points should be sampled by evenly choosing z and phi.

    Args:
        r (float): radius
        rstep (float, optional): radius step.
        azstep (float, optional): azimuthal angle step.
    """
    points = []
    for z in np.arange(-r, r, rstep):
        rz = np.sqrt(r**2 - z**2)
        for az in np.arange(0.0, 2.0*np.pi, azstep):
            # azimuthal angle or rotation around z-axis
            x = rz * np.cos(az)
            y = rz * np.sin(az)
            points.append((x,y,z))

    return np.unique(np.array(points), axis=0)


def rotation_matrix_from_axis_angle(axis:np.ndarray, angle:float) -> np.ndarray:
    """Calculates a rotation matrix from a rotation axis and angle."""
    axis = axis / np.linalg.norm(axis)  # Normalize axis
    c = np.cos(angle)
    s = np.sin(angle)
    x, y, z = axis
    rot_mat = np.array([
        [c + (1 - c) * x**2, (1 - c) * x * y - s * z, (1 - c) * x * z + s * y],
        [(1 - c) * y * x + s * z, c + (1 - c) * y**2, (1 - c) * y * z - s * x],
        [(1 - c) * z * x - s * y, (1 - c) * z * y + s * x, c + (1 - c) * z**2]
    ])

    return rot_mat


def rotate_z_to_point(P:np.ndarray) -> np.ndarray:
    """Calculates rotation matrix to rotate z-axis to a given point on unit sphere."""
    z_axis = np.array([0, 0, 1])
    rotation_axis = np.cross(z_axis, P)
    if np.allclose(rotation_axis, 0):
        if np.allclose(z_axis, P):
            return np.eye(3)
        rotation_axis = np.array([0, 1, 0])
    angle = np.arccos(np.dot(z_axis, P))

    return rotation_matrix_from_axis_angle(rotation_axis, angle)



def uniform_rotation_matrices(rstep:float, azstep:float) -> List[np.ndarray]:
    """Generate evenly distributed rotation matrices.

    In Fast Random Rotation Matrices (James Avro, 1992), 
    a method for uniform random 3D rotation matrices is outlined, the main steps being:

        1. A random rotation about the z axis
        2. Rotate the (unmoved) z unit vector to a random position on the unit sphere
    
    https://www.blopig.com/blog/2021/08/uniformly-sampled-3d-rotation-matrices/

    Returns:
        List[np.ndarray]: a list of rotation matrices
    """
    unit_sphere_points = spherical_points(1.0, rstep, azstep)
    
    U = []
    for az in np.arange(0, 2.*np.pi, azstep):
        R1 = np.array([
            [  np.cos(az), np.sin(az), 0.0 ],
            [ -np.sin(az), np.cos(az), 0.0 ],
            [   0.0      ,  0.0      , 1.0 ]])
        for P in unit_sphere_points:
            R2 = rotate_z_to_point(P)
            U.append(R2 * R1)
    V = np.array(U)
    return np.unique(V, axis=0)


def set_apart(fixed:np.ndarray, 
              movable:np.ndarray,
              r:float = 20.0, 
              steric_clash_cutoff:float=3.0,
              fixed_subidx : Optional[List[int]] = None,
              movable_subidx : Optional[List[int]] = None,
              ) -> Tuple[np.ndarray, float, int]:
    """Find a translation for `movable` coordinates that satisfies the closest distance within desired range.

    Args:
        fixed (np.ndarray): coordinates to be fixed.
        movable (np.ndarry): coordinates to be moved.
        r (float, optional): pocket-to-pocket distance. Defaults to 20.0.

    Returns:
        (np.ndarray,float,int,int): (trans, closest_distance, num_violations, num_trials)
    """
    kd = KDTree(fixed)
    points_on_sphere = spherical_points(r, rstep=2.0, azstep=np.pi/3.0)
    closest_distance = []
    steric_clash = []
    for trans in points_on_sphere:
        coords = movable + trans
        distances, indices = kd.query(coords)
        closest_distance.append(np.min(distances))
        steric_clash = np.sum(distances < steric_clash_cutoff)
    i = np.argmax(closest_distance)
    return (points_on_sphere[i], closest_distance[i], steric_clash[i])


def set_closest_distance(fixed:np.ndarray,
                         movable:np.ndarray, 
                         min_dist:float=5.0,
                         max_dist:float=25.0,
                         max_trials:int=1000) -> Tuple[np.ndarray, float, int, int]:
    """Find a translation for `movable` coordinates that satisfies the closest distance within desired range.

    Args:
        fixed (np.ndarray): coordinates to be fixed.
        movable (np.ndarry): coordinates to be moved.
        min_dist (float, optional): minimum distance allowed. Defaults to 5.0.
        max_dist (float, optional): maximum translation. Defaults to 25.0.
        max_trials (int, optional): maximum number of trials. Defaults to 1000.

    Returns:
        (np.ndarray,float,int,int): (trans, closest_distance, num_violations, num_trials)
    """
    kd = KDTree(fixed)
    distances, indices = kd.query(movable)
    closest_dist = np.min(distances)
    num_violations = np.sum(distances < min_dist) + np.sum(distances > max_dist) # array of True(1) or False(0)
    num_trials = 0
    trans = np.array([0.0, 0.0, 0.0])
    while num_violations > 0 and num_trials < max_trials:
        num_trials += 1
        max_move = abs(0.5*(max_dist + min_dist) - closest_dist)
        trans = max_move * (2.0*np.random.rand(3) -1.0) 
        # numpy.random.rand() generates [0,1) so 2x-1 -> [-1,+1)
        coords = movable + trans
        distances, indices = kd.query(coords)
        num_violations = np.sum(distances < min_dist) + np.sum(distances > max_dist) # array of True(1) or False(0)
        closest_dist = np.min(distances)
    return (trans, closest_dist, num_violations, num_trials)