import numpy as np
import euchar.utils as u 
from scipy.spatial import Delaunay, distance

#=================================================
# Delaunay
#=================================================

def delaunay_ed_tri(points):
    if points.shape[1] != 2:
        raise ValueError("points have invalid number of dimensions")

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
    if points.shape[1] !=3:
        raise ValueError("points have invalid number of dimensions")
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

def dimension(simplex):
    return len(simplex) - np.sum(simplex == -1) - 1

#=================================================

def fix_indices_param(indices, param, simplices):
    for i in range(1, len(indices)):
        p0 = param[i-1]
        p1 = param[i]
        
        if p0 == p1:
            # only need to switch two consecutive simplices
            # alpha par alows only one edge of same parameter
            # as a triangle
            dim_p0 = dimension(simplices[i-1])
            dim_p1 = dimension(simplices[i])
            if dim_p0 > dim_p1:
                higher_dim_simpl = np.copy(simplices[i-1])
                simplices[i-1] = np.copy(simplices[i])
                simplices[i] = higher_dim_simpl
    return simplices

#=================================================

def filter_simplices_on(vertices, simplices, parametrization):

    new_sim = []
    new_par = []

    for k, simpl in enumerate(simplices):
        nope = False
        dim = len(np.flatnonzero(simpl + 1))
        for i in range(dim):
            if simpl[i] not in vertices:
                nope = True
        if nope == False:
            new_sim.append(list(simpl))
            new_par.append(parametrization[k])

    return np.array(new_sim), np.array(new_par)

#=================================================
# Filtrations
#=================================================

def alpha_filtration_2d(points):
    """Compute sorted simplices and sorted parametrization of 2D 
point cloud.
    
    Parameters
    ----------
    points
        array of points in 2D
    
    Returns
    -------
    sorted_simplices
        array of integers of shape (N,3) using -1 as placeholder
        Ex. [1,-1,-1] is vertex 1 and [3,5,-1] is edge (3,5).
    sorted_parametrization
        array of floats corresponding to the filtration parameter of
        each simplex.
    """

    # Make vertices
    vertices = np.arange(len(points))
    edges, triangles = delaunay_ed_tri(points)
    
    # Pamatrise simplices by max(length(edge)) = Rips
    par_vertices = np.zeros(len(points))
    par_edges = np.array([distance.euclidean(*points[ed])
                          for ed in edges])
    par_triangles = np.array([u.parameter_triangle(*points[tri])
                              for tri in triangles])
    
    # Unsorted simplices using arrays with same sizes.
    vertices = np.array([[v, -1, -1] for v in vertices])
    edges = np.array([[ed[0], ed[1], -1] for ed in edges])
    unsorted_simplices = np.vstack([vertices, edges, triangles])
    unsorted_param = np.hstack([par_edges, par_triangles])

    # Make sorted indices (knowing all parameters of vertices are 0)
    sorted_indices_edges_triangles = np.argsort(unsorted_param)
    sorted_indices_edges_triangles += len(points)
    sorted_indices = np.hstack([np.arange(len(vertices)),
                                sorted_indices_edges_triangles])

    # Sort both simplices and their param with sorted_indices
    sorted_simplices = unsorted_simplices[sorted_indices]
    unsorted_param = np.hstack([np.zeros(len(vertices)),
                                unsorted_param])
    sorted_param = unsorted_param[sorted_indices]

    # Fix order edges and triangles
    sorted_simplices = fix_indices_param(sorted_indices,
                                         sorted_param,
                                         sorted_simplices)

    return sorted_simplices, sorted_param

#=================================================

def alpha_filtration_3d(points):
    """
    Return sorted parametrization and sorted simplices.
    
    Parameters
    ----------
    points
        array of points in 3D
    
    Returns
    -------
    sorted_param
        array of floats corresponding to the filtration parameter of
        each simplex.
    sorted_simplices
        array of integers of shape (N,4) using -1 as placeholder
        Ex. [1,-1,-1, -1] is vertex 1 and [3,5,-1, -1] is edge (3,5).
    """
    
    vertices = np.arange(len(points))
    edges, triangles, tetrahedra = delaunay_ed_tri_tetra(points)

    # Pamatrise simplices by max(length(edge)) = Rips
    par_vertices = np.zeros(len(points))
    par_edges = np.array([distance.euclidean(*points[ed])
                          for ed in edges])
    par_triangles = np.array([u.parameter_triangle(*points[tri])
                              for tri in triangles])
    
    par_tetrahedra = np.array([u.parameter_tetrahedron(*points[tetra])
                           for tetra in tetrahedra])

    # Unsorted simplices using arrays with same sizes. 
    vertices = np.array([[v, -1, -1, -1] for v in vertices])
    edges = np.array([[ed[0], ed[1], -1, -1] for ed in edges])
    triangles = np.array([[tri[0], tri[1], tri[2], -1] for tri in triangles])
    
    unsorted_simplices = np.vstack([vertices, edges, triangles, tetrahedra])
    unsorted_param = np.hstack([par_edges, par_triangles, par_tetrahedra])

    # Make sorted indices (knowing all parameters of vertices are 0)
    sorted_indices_edges_triangles_tetrahedra = np.argsort(unsorted_param)
    sorted_indices_edges_triangles_tetrahedra += len(points)
    sorted_indices = np.hstack([np.arange(len(vertices)), sorted_indices_edges_triangles_tetrahedra])

    # Sort both simplices and their parametrization with sorted_indices
    sorted_simplices = unsorted_simplices[sorted_indices]
    unsorted_param = np.hstack([np.zeros(len(vertices)), unsorted_param])
    sorted_param = unsorted_param[sorted_indices]

    # fix order edges and triangles
    sorted_simplices = fix_indices_param(sorted_indices, sorted_param, sorted_simplices)

    return sorted_simplices, sorted_param

#=================================================

def density_filtration(points, simplices, n_neighbors):
    """Parametrization of simplices based on local density at vertices."""

    density = u.estimate_density(points, n_neighbors)

    parametrization = []
    
    for sim in simplices:
        par = max([density[vertex] for vertex in sim if vertex != -1])
        parametrization.append(par)

    return np.array(parametrization)
