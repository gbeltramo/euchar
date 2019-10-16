import numpy as np
from scipy.spatial.distance import euclidean
from sklearn.neighbors import NearestNeighbors
from euchar.cppbinding.utils import vector_of_euler_changes_2d, vector_of_euler_changes_3d

#=================================================

def estimate_inverse_density(points, n):
    """
    Estimates the inverse of the density at each `p` in `points` as the
    root mean square of the distances from `p` to its `n` nearest 
    neighbors.

    Parameters
    ----------
    points
        np.ndarray of shape (N, d)
    n
        int, number of nearest neighbors used to estimate the inverse of
        the density.
    
    Returns
    -------
    param_vertices
        1-D np.ndarray parametrizing vertices by inverse density.

    """
    nbrs = NearestNeighbors(n+1, algorithm="ball_tree")
    dist = nbrs.fit(points).kneighbors(points)[0][:,1:]

    densities = np.empty(len(points))
    for i in range(len(points)):
        densities[i] = np.sqrt(np.sum((dist[i]**2) / len(dist[i])))
    return densities

#=================================================

def vector_all_euler_changes_in_2D_images():
    """
    Returns a vector of 256 integers, containing all the possible Euler
    characteristic changes produced by the insertion of a central pixel
    in a 3x3 binary neighbourhood.
    """
    
    list_2d_changes = vector_of_euler_changes_2d()
    return np.array(list_2d_changes)

#=================================================

def vector_all_euler_changes_in_3D_images(start=0, end=67_108_864):
    """
    Returns a vector of 67,108,864 integers, containing all the
    possible Euler characteristic changes produced by the insertion of
    a central voxel in a 3x3x3 binary neighbourhood.
    
    Notes
    -----
    It is recommended to produce this vector once and to write it 
    to a file. For example with numpy

    >>> vector_changes = vector_all_euler_changes_in_3d_images()
    >>> np.save("/path/to/dir/changes_3d.npy", vector_changes)
    
    The next time the vector can be loaded with

    >>> vector_changes = np.load("/path/to/dir/changes_3d.npy")
    """
    
    print("\n-----\nIt is recommended to save the result to a file, and to load it when needed.\n-----\n")

    list_3d_changes = vector_of_euler_changes_3d(start, end)
    return np.array(list_3d_changes)

#=================================================

def magnitude(v):
    """Euclidean norm of the vector v."""
    return np.sqrt(v.dot(v))

#=================================================

def circumradius_triangle(A, B, C):
    """
    Compute circumradius of 3 points in 2 or 3 dimensions.
    
    Parameters
    ----------
    A, B, C
        np.arrays of length equal to the Euclidean dimension.

    Return
    ------
    circumradius
        float, circumradius of the 3 points.
    """
    a = np.array(A) - np.array(C)
    b = np.array(B) - np.array(C)
    cross = np.cross(a,b)
    return magnitude(a) * magnitude(b) * magnitude(a-b) \
           / 2 / magnitude(cross)

#=================================================

def circumcenter_triangle(A, B, C):
    """
    Compute circumcenter of 3 points in 2 or 3 dimensions.

    Parameters
    ----------
    A, B, C
        np.arrays of length equal to the Euclidean dimension.

    Return
    ------
    circumcenter
        np.array, the circumcenter of the three points
    """
    a = np.array(A) - np.array(C)
    b = np.array(B) - np.array(C)
    aXb = np.cross(a,b)
    numer1 = magnitude(a)**2 * b - magnitude(b)**2 * a
    numer = np.dot(numer1, b)*a - np.dot(numer1, a)*b
    denom = 2 * (magnitude(aXb)**2)
    return numer / denom + C

#=================================================

def vol_tetra(A, B, C, D):
    """
    Volume of the tetrahedron whose four vertices are the given four 
    points A, B, C, D.

    Parameters
    ----------
    A, B, C, D
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    volume
        float, the volume of the tetrahedron
    """
    
    a = A - D
    b = B - D
    c = C - D
    volume = abs(np.dot(a, np.cross(b,c))) / 6
    return volume

#=================================================

def circumradius_tetra(A, B, C, D):
    """
    Compute circumradius of 4 points in 3 dimensions.

    Parameters
    ----------
    A, B, C, D
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    circumradius
        np.array, the circumradius of the four points
    """
    
    e12 = magnitude(A - B) # length edge AB
    o12 = magnitude(C - D) # length edge opposite to AB
    prod12 = e12 * o12
    
    e13 = magnitude(A - C) # same
    o13 = magnitude(B - D)
    prod13 = e13 * o13
    
    e14 = magnitude(A - D) # same
    o14 = magnitude(B - C)
    prod14 = e14 * o14
    
    numer = np.sqrt((prod12 + prod13 + prod14)*
                    (prod12 + prod13 - prod14)*
                    (prod12 - prod13 + prod14)*
                    (-prod12 + prod13 + prod14))
    denom = 24 * vol_tetra(A, B, C, D)
    circumradius = numer / denom
    return circumradius

#=================================================

def circumcenter_tetra(A, B, C, D):
    """
    Compute circumcenter of 4 points in 3 dimensions.

    Parameters
    ----------
    A, B, C, D
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    circumcenter
        np.array, the circumcenter of the four points
    """
    
    
    colA1, colA2, colA3 = A - D, B - D, C - D
    matrix = np.vstack([colA1, colA2, colA3])
    inverse_matrix = np.linalg.inv(matrix)
    
    B = np.array([A.dot(A) - D.dot(D), B.dot(B) - D.dot(D),
                  C.dot(C) - D.dot(D)]) / 2
    
    circumcenter = [np.dot(row, B) for row in inverse_matrix]
    return circumcenter

#=================================================

def is_acute(lengths):
    """Check if a triangle is acute with Pythgoras theorem."""
    l = sorted(lengths)
    return l[0]*l[0] + l[1]*l[1] >= l[2]*l[2]

#=================================================

def parameter_triangle(A, B, C):
    """
    Compute the Alpha filtration parameter of a triangle of vertices 
    A, B, C.
    
    Parameters
    ----------
    A, B, C
        np.arrays, the three vertices of the triangle

    Return
    ------
    parameter
        float, Alpha filtration parameter of the triangle. This
        is twice the circumradius of the triangle if the triangle is
        acute, or it is the greatest value among the length of edges 
        of the triangle.

    """
    lengths = [euclidean(A, B), euclidean(B, C), euclidean(A, C)]
    acute = is_acute(lengths)
    
    if acute:
        parameter = 2 * circumradius_triangle(A, B, C)
    else:
        parameter = np.array([euclidean(A, B),
                              euclidean(B, C),
                              euclidean(A, C)]).max()
    
    return parameter


#=================================================

def center_triangle(A, B, C):
    lengths = [euclidean(A, B), euclidean(B, C), euclidean(A, C)]
    acute = is_acute(lengths)
    
    if acute:
        C = circumcenter_triangle(A, B, C)
    else:
        ind_par = np.argmax(lengths)
        centers = np.array([(A + B) / 2, 
                            (B + C) / 2,
                            (A + C) / 2])
        C = centers[ind_par]
    
    return C

#=================================================

def is_on_correct_side(A, B, C, test_point, point):
    cr = np.cross(B-A, C-A)
    dot_test = np.dot(cr, test_point-A)
    dot_point = np.dot(cr, point-A)
    
    if dot_test >= 0 and dot_point >= 0:
        return True
    elif dot_test <= 0 and dot_point <= 0:
        return True
    else:
        return False

#=================================================

def circumcenter_is_inside(circum, A, B, C, D):
    if not is_on_correct_side(A, B, C, test_point=D, point=circum):
        return False
    
    if not is_on_correct_side(A, B, D, test_point=C, point=circum):
        return False
    
    if not is_on_correct_side(A, C, D, test_point=B, point=circum):
        return False
    
    if not is_on_correct_side(C, B, D, test_point=A, point=circum):
        return False
    
    return True

#=================================================
    
def parameter_tetrahedron(A, B, C, D):
    """
    Compute the Alpha filtration parameter of a tetrahedron of vertices
    A, B, C, D.
    
    Parameters
    ----------
    A, B, C, D
        np.arrays, the four vertices of the tetrahedron.

    Return
    ------
    parameter
        float, Alpha filtration parameter of the tetrahedron.
    
    """

    circumcenter = circumcenter_tetra(A, B, C, D)

    if circumcenter_is_inside(circumcenter, A, B, C, D):
        return 2 * circumradius_tetra(A, B, C, D)
    else:
        par_ABC = parameter_triangle(A, B, C)
        par_ABD = parameter_triangle(A, B, D)
        par_ACD = parameter_triangle(A, C, D)
        par_BCD = parameter_triangle(B, C, D)
        return max(par_ABC, par_ABD, par_ACD, par_BCD)

#=================================================

def center_tetrahedron(A, B, C, D):
    circum = circumcenter_tetra(A, B, C, D)

    if circumcenter_is_inside(circum, A, B, C, D):
        return circum
    else:
        C1 = center_triangle(A, B, C)
        C2 = center_triangle(A, B, D)
        C3 = center_triangle(A, C, D)
        C4 = center_triangle(B, C, D)
    
        centers = [C1, C2, C3, C4]
    
        par_triABC = parameter_triangle(A, B, C)
        par_triABD = parameter_triangle(A, B, D)
        par_triACD = parameter_triangle(A, C, D)
        par_triBCD = parameter_triangle(B, C, D)
    
        parameters = [par_triABC, par_triABD, par_triACD, par_triBCD]
        ind = np.argmax(np.array(parameters))
        return centers[ind]
    
#=================================================

def simplices_to_dimensions(simplices):
    """
    Transforms an array of simplices to an array of their dimensions.

    Example
    -------
    Given the array of simplices
    
    >>> [[1 -1 -1], [3 4 -1], [10 12 15]]

    the output will be 

    >>> [0 1 2]

    """
    def dim_simplex(sim):
        cnt = 0
        for el in sim:
            if el != -1:
                cnt += 1
            else:
                break
        return cnt - 1

    return np.apply_along_axis(dim_simplex, 1, simplices)
