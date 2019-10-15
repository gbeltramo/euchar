import numpy as np
import euchar.utils as u 
from scipy.spatial import Delaunay
from scipy.spatial.distance import euclidean

#=================================================

def delaunay_ed_tri(points):
    """Delaunay triangulation simplices of two dimensional points."""
    if points.shape[1] != 2:
        raise ValueError("Points need to be a two dimensional np.ndarray.")

    tri = Delaunay(points)
    arr_triangles = np.array([np.sort(unsort_tri)
                              for unsort_tri in tri.simplices])

    edges = set()
    for t in arr_triangles:
        i1, j1 = t[0], t[1]
        i2, j2 = t[0], t[2]
        i3, j3 = t[1], t[2]
        edges.add((i1, j1))
        edges.add((i2, j2))
        edges.add((i3, j3))

    return np.array(list(edges)), arr_triangles

#=================================================

def delaunay_ed_tri_tetra(points):
    """Delauanay triangulation simplices of three dimensional points."""
    if points.shape[1] !=3:
        raise ValueError("Points need to be a three dimensional np.ndarray.")
    tetra = Delaunay(points)
    arr_tetrahedra = np.array([np.sort(unsort_tetra)
                               for unsort_tetra in tetra.simplices])
    edges = set()
    triangles = set()
    for tetra in arr_tetrahedra:
        v1, v2, v3, v4 = tetra
        edges.add((v1, v2))
        edges.add((v1, v3))
        edges.add((v2, v3))
        triangles.add((v1, v2, v3))
        triangles.add((v1, v2, v4))
        triangles.add((v1, v3, v4))
        triangles.add((v2, v3, v4))

    return (np.array(list(edges)), np.array(list(triangles)),
                    arr_tetrahedra)

#=================================================

def alpha_filtration_2D(points):
    """
    Simplices and corresponding parametrization of 2D point cloud.
    
    Parameters
    ----------
    points
        np.ndarray of shape (N, 2)
    
    Returns
    -------
    simplices
        array of integers of shape (N,3) using -1 as placeholder
        Ex. [1,-1,-1] is vertex 1 and [3,5,-1] is edge (3,5).
    parametrization
        array of floats corresponding to the Alpha filtration parameter
        of each simplex.
    
    """

    # Get vertices, and Delaunay simplices
    vertices = np.arange(len(points))
    edges, triangles = delaunay_ed_tri(points)

    # Make simplices, with -1 as a placeholder
    vertices = np.array([[v, -1, -1] for v in vertices])
    edges = np.array([[ed[0], ed[1], -1] for ed in edges])
    simplices = np.vstack([vertices, edges, triangles])
    
    # Parametrization
    par_vertices = np.zeros(len(points))
    par_edges = np.array([euclidean(points[ed[0]],
                                    points[ed[1]])
                          for ed in edges])
    par_triangles = np.array([u.parameter_triangle(points[tri[0]],
                                                   points[tri[1]],
                                                   points[tri[2]])
                              for tri in triangles])
    parametrization = np.hstack([par_vertices, par_edges,
                                 par_triangles])

    return simplices, parametrization

#=================================================

def alpha_filtration_3D(points):
    """
    Simplices and corresponding parametrization of 3D point cloud.
    
    Parameters
    ----------
    points
        np.ndarray of shape (N, 3)
    
    Returns
    -------
    simplices
        array of integers of shape (N, 4) using -1 as placeholder
        Ex. [1,-1,-1,-1] is vertex 1 and [3,5,-1,-1] is edge (3,5).
    parametrization
        array of floats corresponding to the Alpha filtration parameter
        of each simplex.
    
    """
    # Get vertices, and Delaunay simplices
    vertices = np.arange(len(points))
    edges, triangles, tetrahedra = delaunay_ed_tri_tetra(points)

    # Make simplices, with -1 as a placeholder
    vertices = np.array([[v, -1, -1, -1] for v in vertices])
    edges = np.array([[ed[0], ed[1], -1, -1] for ed in edges])
    triangles = np.array([[tri[0], tri[1], tri[2], -1]
                          for tri in triangles])
    simplices = np.vstack([vertices, edges, triangles, tetrahedra])
    
    # Parametrization
    par_vertices = np.zeros(len(points))
    par_edges = np.array([euclidean(points[ed[0]],
                                    points[ed[1]])
                          for ed in edges])
    par_triangles = np.array([u.parameter_triangle(points[tri[0]],
                                                   points[tri[1]],
                                                   points[tri[2]])
                              for tri in triangles])
    par_tetrahedra = np.array([u.parameter_tetrahedron(*points[tetra])
                               for tetra in tetrahedra])
    parametrization = np.hstack([par_vertices, par_edges,
                                 par_triangles, par_tetrahedra])

    return simplices, parametrization

#=================================================

def inverse_density_filtration(points, simplices, n_neighbors):
    """
    Parametrization of simplices based on inverse of local density
    at vertices. The inverse of the density at each `p` in `points`
    is estimated as the root mean square of the distances from `p` to
    its `n` nearest neighbors.
    
    Parameters
    ----------
    points
       np.ndarray of shape (N, d)
    simplices
        np.ndarray of shape (M, d+1)
    n_neighbors
        int, number of nearest neighbors of a point `p` in `points`
        used to estimate the local density around `p`.
    
    Return
    ------
    density_parametrization
        1-D np.ndarray, parametrization of `simplices`. The value of a 
        simplex is the maximum of any of its vertices.
    """

    density = u.estimate_inverse_density(points, n_neighbors)

    parametrization = []
    
    for sim in simplices:
        par = max([density[vertex] for vertex in sim if vertex != -1])
        parametrization.append(par)

    return np.array(parametrization)
