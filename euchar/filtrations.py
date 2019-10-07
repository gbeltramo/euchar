import numpy as np
from scipy.spatial import Delaunay, distance

#=================================================

def delaunay_edges_and_triangles(points):

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

    return np.array(list(edges)), np.array(list(triangles)), arr_tetrahedra

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

def smallest_circle2d(a, b, c):
    c12x = (a[0]+b[0])/2
    c12y = (a[1]+b[1])/2
    c12 = [c12x, c12y]
    r12 = np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)/2

    c13x = (a[0]+c[0])/2
    c13y = (a[1]+c[1])/2
    c13 = [c13x, c13y]
    r13 = np.sqrt((a[0] - c[0])**2 + (a[1] - c[1])**2)/2

    c23x = (b[0]+c[0])/2
    c23y = (b[1]+c[1])/2
    c23 = [c23x, c23y]
    r23 = np.sqrt((b[0] - c[0])**2 + (b[1] - c[1])**2)/2

    dist12_3 = distance.euclidean(c12, c)
    dist13_2 = distance.euclidean(c13, b)
    dist23_1 = distance.euclidean(c23, a)
    list_dist = [dist12_3, dist13_2, dist23_1]

    list_c = np.array([c12, c13, c23])
    list_r = np.array([r12, r13, r23])
    list_ok = [(dist12_3 < r12), (dist13_2 < r13), (dist23_1 < r23)]

    if sum(list_ok) != 0:
        #print("list_c", list_c)
        #print("list_r", list_r)
        #print("2list_r", 2*list_r)
        #print("list_dist", list_dist)
        #print("list_ok", list_ok)
        ind = np.argmin(list_r[list_ok])
        return list_c[list_ok][ind], list_r[list_ok][ind]

    # translate of A  : https://en.wikipedia.org/wiki/Circumscribed_circle
    bx = b[0] - a[0]
    by = b[1] - a[1]
    cx = c[0] - a[0]
    cy = c[1] - a[1]

    D = 2*(bx*cy - by*cx)
    ux = 1/D *( cy*(bx**2 + by**2) - by*(cx**2 + cy**2))
    uy = 1/D *( bx*(cx**2 + cy**2) - cx*(bx**2 + by**2))
    rrr = np.sqrt(ux**2 + uy**2)
    c123x = ux+a[0]
    c123y = uy+a[1]

    return [c123x, c123y], rrr


#=================================================

def delaunay_filtration_2d(points, print_info=False):
    """
    Return tuple of arrays of sorted simplices and sorted parametrization
    
    Parameters
    ----------
    pts
        array of points in 2D
    
    Returns
    -------
    sorted_simplices
        array of integers of shape (N,3) using -1 as placeholder
        Ex. [1,-1,-1] is vertex 1 and [3,5,-1] is edge (3,5).
    sorted_parametrization
        array of floats corresponding to the filtration parameter of
        each simplex.
        
    Notes
    -----
    Implementation uses: scipy.spatial.Delaunay, scipy.spatial.distance,
    tda.delaunay.smallest_circle2d. 
    After finding Delaunay edges and triangles we use np.argsort to 
    sort both simplices and their parametrization.
    """

    # Make vertices
    arr_vertices = np.arange(len(points))
    arr_edges, arr_triangles = delaunay_edges_and_triangles(points)

    if print_info:
        print("edges and triangles done", end=" - ")
        
    # Pamatrise simplices by max(length(edge)) = Rips
    par_vertices = np.zeros(len(points))
    par_edges = np.array([distance.euclidean(points[ed[0]], points[ed[1]]) for ed in arr_edges])
    dict_edges = {tuple(ed): distance.euclidean(points[ed[0]], points[ed[1]]) for ed in arr_edges}
    par_triangles = np.array([max(dict_edges[(i1,i2)],
                                  dict_edges[(i1,i3)],
                                  dict_edges[(i2,i3)],
                                  2*smallest_circle2d(points[i1], points[i2], points[i3])[1])
                              for i1, i2, i3 in arr_triangles])
    if print_info:
        print("paramaetrizations done", end=" - ")
        
    # Unsorted simplices using arrays with same sizes. -1 is a placeholder
    arr_vertices = np.array([[v, -1, -1] for v in arr_vertices])
    arr_edges = np.array([[ed[0], ed[1], -1] for ed in arr_edges])
    unsorted_simplices = np.vstack([arr_vertices, arr_edges, arr_triangles])
    unsorted_parametrization = np.hstack([par_edges, par_triangles])

    # Make sorted indices (knowing all parameters of vertices are 0)
    sorted_indices_edges_triangles = np.argsort(unsorted_parametrization)
    sorted_indices_edges_triangles += len(points)
    sorted_indices = np.hstack([np.arange(len(arr_vertices)), sorted_indices_edges_triangles])

    # Sort both simplices and their parametrization with sorted_indices
    sorted_simplices = unsorted_simplices[sorted_indices]
    unsorted_parametrization = np.hstack([np.zeros(len(arr_vertices)), unsorted_parametrization])
    sorted_parametrization = unsorted_parametrization[sorted_indices]

    # fix order edges and triangles
    sorted_simplices = fix_indices_param(sorted_indices, sorted_parametrization, sorted_simplices)

    return sorted_simplices, sorted_parametrization

