import numpy as np
from scipy.spatial.distance import euclidean
from sklearn.neighbors import NearestNeighbors
from euchar.cppbinding.utils import vector_of_euler_changes_2d, vector_of_euler_changes_3d

#=================================================

def estimate_density(points, n_neighbors, algorithm="ball_tree"):
    """
    Estimates density at each point in `points`.

    Parameters
    ----------
    points
        2d array of shape (N, d)
    n_neighbors
        number of closest points to consider when estimating the density
    algorithm
        string specifying which algorithm the NearestNeighbors class
        from sklearn.neighbors should use

    Returns
    -------
    param_vertices
        1d array parametrizing vertices by density. The order is the one of
        indices of vertices (0, 1, 2, ...), so parametrized values are unsorted

    Notes
    -----
    This requires scikit-learn to be installed in the current 
    environment.
    """
    nbrs = NearestNeighbors(n_neighbors+1, algorithm)
    dist = nbrs.fit(points).kneighbors(points)[0][:,1:]

    densities = np.empty(len(points))
    for i in range(len(points)):
        densities[i] = np.sqrt(np.sum((dist[i]**2) / len(dist[i])))
    return densities

#=================================================

def vector_all_euler_changes_in_2D_images():
    """Returns a vector of 256 integers, containing all the possible Euler characteristic changes produced by the insertion of a central pixel in a 3x3 binary neighbourhood.
    """
    
    list_2d_changes = vector_of_euler_changes_2d()
    return np.array(list_2d_changes)

#=================================================

def vector_all_euler_changes_in_3D_images():
    """Returns a vector of 67,108,864 integers, containing all the possible Euler characteristic changes produced by the insertion of a central pixel in a 3x3x3 binary neighbourhood.
    
    Notes
    -----
    It is recommended to produce this vector once and to write it 
    to a file. For example with numpy

    >>> vector_changes = vector_all_euler_changes_in_3d_images()
    >>> np.save("/path/to/dir/changes_3d.npy", vector_changes)
    
    The next time the vector can be loaded with

    >>> vector_changes = np.load("/path/to/dir/changes_3d.npy")
    """
    
    print("This may take a while...")
    print("It is recommended to save the result to a file,\n and to load it when needed.")

    list_3d_changes = vector_of_euler_changes_3d(0, 67_108_864) # 2**26
    return np.array(list_3d_changes)

#=================================================

def magnitude(v):
    """Euclidean norm of the vector v."""
    return np.sqrt(v.dot(v))

#=================================================

def circumradius_triangle(A, B, C):
    """Compute circumradius of 3 points in 2 or 3 dimensions.
    
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
    """Compute circumcentre of 3 points in 2 or 3 dimensions.

    Parameters
    ----------
    A, B, C
        np.arrays of length equal to the Euclidean dimension.

    Return
    ------
    circumcentre
        np.array, the circumcentre of the three points
    """
    a = np.array(A) - np.array(C)
    b = np.array(B) - np.array(C)
    aXb = np.cross(a,b)
    numer1 = magnitude(a)**2 * b - magnitude(b)**2 * a
    numer = np.dot(numer1, b)*a - np.dot(numer1, a)*b
    denom = 2 * (magnitude(aXb)**2)
    return numer / denom + C

#=================================================

def vol_tetra(p1, p2, p3, p4):
    """Volume of the tetrahedron whose four vertices are the given four points p1, p2, p3, p4.

    Parameters
    ----------
    p1, p2, p3, p4
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    volume
        float, the volume of the tetrahedron
    """
    
    a = p1 - p4
    b = p2 - p4
    c = p3 - p4
    volume = abs(np.dot(a, np.cross(b,c))) / 6
    return volume

#=================================================

def circumradius_tetra(p1, p2, p3, p4):
    """Compute circumradius of 4 points in 3 dimensions.

    Parameters
    ----------
    p1, p2, p3, p4
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    circumradius
        np.array, the circumradius of the four points
    """
    
    e12 = magnitude(p1 - p2) # length edge p1p2
    o12 = magnitude(p3 - p4) # length edge opposite to p1p2
    prod12 = e12 * o12
    
    e13 = magnitude(p1 - p3) # same
    o13 = magnitude(p2 - p4)
    prod13 = e13 * o13
    
    e14 = magnitude(p1 - p4) # same
    o14 = magnitude(p2 - p3)
    prod14 = e14 * o14
    
    numer = np.sqrt((prod12 + prod13 + prod14)*
                    (prod12 + prod13 - prod14)*
                    (prod12 - prod13 + prod14)*
                    (-prod12 + prod13 + prod14))
    denom = 24 * vol_tetra(p1, p2, p3, p4)
    circumradius = numer / denom
    return circumradius

#=================================================

def circumcenter_tetra(p1, p2, p3, p4):
    """Compute circumcenter of 4 points in 3 dimensions.

    Parameters
    ----------
    p1, p2, p3, p4
        np.arrays, the four vertices of the tetrahedron

    Return
    ------
    circumcenter
        np.array, the circumcentre of the four points
    """
    
    
    colA1, colA2, colA3 = p1 - p4, p2 - p4, p3 - p4
    A = np.vstack([colA1, colA2, colA3])
    Ainv = np.linalg.inv(A)
    
    B = np.array([p1.dot(p1) - p4.dot(p4), p2.dot(p2) - p4.dot(p4),
                  p3.dot(p3) - p4.dot(p4)]) / 2
    
    circumcentre = [np.dot(row, B) for row in Ainv]
    return circumcentre

#=================================================

def parameter_triangle(A, B, C):
    """Compute the Alpha filtration parameter of a triangle of vertices A, B, C.
    
    Parameters
    ----------
    A, B, C
        np.arrays, the three vertices of the triangle

    Return
    ------
    parameter
        float, Alpha filtration parameter of the triangle. This
        is twice the circumradius of the triangle if this is
        greater than of all the three triangle edges lengths, or
        it is the greatest value among the length of edges of the
        triangle.

    """
    
    AB = euclidean(A, B)
    AC = euclidean(A, C)
    BC = euclidean(B, C)
    twice_circ_radius_ABC = 2 * circumradius_triangle(A, B, C)
    return max(AB, AC, BC, twice_circ_radius_ABC)

#=================================================
    
def parameter_tetrahedron(p1, p2, p3, p4):
    """Compute the Alpha filtration parameter of a tetrahedron of vertices p1, p2, p3, p4.
    
    Parameters
    ----------
    p1, p2, p3, p4
        np.arrays, the four vertices of the tetrahedron.

    Return
    ------
    parameter
        float, Alpha filtration parameter of the tetrahedron. This
        is the greatest value among twice the circumradius of the 
        tetrahedron and twice the circumradiuses of the four 
        triangles forming the faces of the tetrahedron.

    """
    
    twice_circ_radius_tetra = 2 * circumradius_tetra(p1, p2, p3, p4)
    twice_circ_radius_tri123 = 2 * circumradius_triangle(p1, p2, p3)
    twice_circ_radius_tri124 = 2 * circumradius_triangle(p1, p2, p4)
    twice_circ_radius_tri134 = 2 * circumradius_triangle(p1, p3, p4)
    twice_circ_radius_tri234 = 2 * circumradius_triangle(p2, p3, p4)
    
    return max(twice_circ_radius_tetra, twice_circ_radius_tri123,
               twice_circ_radius_tri124, twice_circ_radius_tri134,
               twice_circ_radius_tri234)

#=================================================

def simplices_to_dimensions(simplices):

    def dim_simplex(sim):
        cnt = 0
        for el in sim:
            if el != 1:
                cnt += 1
            else:
                break
        return cnt - 1

    return [dim_simplex(item) for item in simplices]

        
