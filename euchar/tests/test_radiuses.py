import numpy as np
from euchar.utils import circumradius_triangle, circumcenter_triangle, circumcenter_tetra, circumcenter_tetra

#=================================================

def test_circumradius_triangle_2d(N=100):
    triangles = np.random.rand(N, 3, 2)
    radiuses = [circumradius_triangle(*tri) for tri in triangles]
    centers = [circumcenter_triangle(*tri) for tri in triangles]

    distances = [np.isclose(distance.euclidean(tri[0], c), r)
                 for tri, c, r in zip(triangles, centers, radiuses)]
    assert np.all(distances), "Alpha filtration radius of 2D triangles does not match distance from circumcentre."

#=================================================

def test_circumradius_tetra(N=100):
    tetrahedra = np.random.rand(N, 4, 3)
    radiuses = [circumradius_tetra(*tetra) for tetra in tetrahedra]
    centers =  [circumcenter_tetra(*tetra) for tetra in tetrahedra]

    distances = [np.isclose(distance.euclidean(tetra[0], c), r)
                 for tetra, c, r in zip(tetrahedra, centers, radiuses)]
    assert np.all(distances), "Alpha filtration radius of 3D tetrahedra does not match distance form circumcentre."
