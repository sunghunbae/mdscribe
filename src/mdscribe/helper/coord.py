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



def avoid_clash(fixed:np.ndarray, movable:np.ndarray, max_move:float=25.0, clash_cutoff:float=5.0, max_trials:int=1000) -> Tuple[int, int, float, np.ndarray]:
    """Find a translation vector for `movable` coordinates to avoid steric clashes.

    Args:
        fixed (np.ndarray): coordinates to be fixed.
        movable (np.ndarry): coordinates to be moved.
        max_move (float, optional): maximum translation. Defaults to 25.0.
        clash_cutoff (float, optional): distance considered as potential clash. Defaults to 5.0.
        max_trials (int, optional): maximum number of trials. Defaults to 1000.

    Returns:
        (int,int,float,np.ndarray): (num_trials, clash_count, closest_dist, trans)
    """
    kd = KDTree(fixed)
    num_trials = 0
    clash_count = 1
    while clash_count > 0 and num_trials < max_trials:
        num_trials += 1
        trans = max_move * (2.0*np.random.rand(3) -1.0) 
        # numpy.random.rand() generates [0,1) so 2x-1 -> [-1,+1)
        coords = movable + trans
        distances, indices = kd.query(coords)
        clash_count = np.sum(distances < clash_cutoff)
        closest_dist = np.min(distances)
    return (num_trials, clash_count, closest_dist, trans)