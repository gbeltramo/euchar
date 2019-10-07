import numpy as np
from sklearn.neighbors import NearestNeighbors

#=================================================

def estimate_density(points, n_neighbors, algorithm="ball_tree"):
    """
    Estimates density at each point in points

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
    """
    nbrs = NearestNeighbors(n_neighbors+1, algorithm)
    dist = nbrs.fit(points).kneighbors(points)[0][:,1:]

    densities = np.empty(len(points))
    for i in range(len(points)):
        densities[i] = np.sqrt(np.sum((dist[i]**2) / len(dist[i])))
    return densities

#=================================================
