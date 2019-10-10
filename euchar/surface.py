"""


    Implementation of the Euler characateristic surfaces algorithms.

    Author: Gabriele Beltramo

"""

import numpy as np
import euchar.utils as u
import euchar.filtrations as f
import euchar.curve as c
import euchar.cppbinding.surface as cppsurface

#=================================================

def images_2D(image1, image2, vector_of_euler_changes_2D=None,
              max_intensity1=255, max_intensity2=255):
    """Euler characteristic curve of a pair of 2D images.
    
    This uses the vector of all possible Euler characteristic changes 
    produced by a pixel insertion in 2D images. This is recomputed 
    every time this function is called if it is not passes as a 
    parameter.

    Parameters
    ----------
    image1, image2
        np.ndarray of integers
    vector_of_euler_changes_2D
        list of integers, precomputed Euler characteristic changes
        produced by a single pixel insertion
    max_intensity1, max_intensity2
        maximum value of any input of the form of image1 and image2
    
    Return
    ------
    euler_char_surface
        np.ndarray of integers

    """

    try:
        err_line1 = "`image` parameter needs to be an np.ndarray\n---"
        assert str(type(image1)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(image2)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a two dimensional np.array\n---"
        assert len(image1.shape) == 2, err_line1
        assert len(image2.shape) == 2, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    if vector_of_euler_changes_2D is None:
        vector_of_euler_changes_2D = u.vector_all_euler_changes_in_2D_images()

    euler_char_surface = cppsurface.images_2d(image1, image2,
                                            vector_of_euler_changes_2D,
                                            max_intensity1,
                                            max_intensity2)

    return np.array(euler_char_surface)

#=================================================

def image_with_function(image, iterations, func, max_intensity=255,
                        **kwargs):
    """Compute the Euler characteristic surface of a give image with the
    bi-filtration given by pixel value and iterative application of
    a given function.

    Parameters
    ----------
    image
        np.ndarray of integers
    iterations
        number of times the function is applied to the image.
    func
        function taking an image at returning an image of the
        same size. The default is cv2.GaussianBlur, so you
        need to pass something like `**{'ksize':(3,3), 'sigmaX':1.0}`
        if you do not set this to something else.
    max_intensity
        int, maximum pixel value.
    **kwargs
        possible arguments of function `func` passed as a parameter.

    Returns
    -------
    euler_char_surface
        np.ndarray of integers

    """
    
    surface = np.empty((max_intensity+1, iterations+1), dtype=int)

    # Copy `image` into `new_image` to avoid aliasing
    new_image = np.copy(image)

    for j in range(iterations+1):
        # We use iterations+1 so to have also the
        # original `image` euler curve
        vector = u.vector_all_euler_changes_in_2D_images()
        ecc = c.image_2d(new_image, vector, max_intensity)
        surface[:,j] = np.array(ecc)

        new_image = func(image, **kwargs)

    return surface


#=================================================

def images_3D(image1, image2, vector_of_euler_changes_3D=None,
              max_intensity1=255, max_intensity2=255):
    """Euler characteristic surface of 3D images.
    
    This uses the vector of all possible Euler characteristic changes 
    produced by a voxel insertion in 3D images. This is must be passed
    as a parameter.
    
    Parameters
    ----------
    image1, image2
        np.ndarray of integers
    vector_of_euler_changes_3D
        list of integers, precomputed Euler characteristic changes
        produced by a single voxel insertion
    max_intensity1, max_intensity2
        maximum value of any input of the form of image1 and image2
    
    Return
    ------
    euler_char_surface
        np.ndarray of integers

    """

    try:
        err_line1 = "`image` parameter needs to be an np.array\n---"
        assert str(type(image1)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(image2)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a three dimensional np.array\n---"
        assert len(image1.shape) == 3, err_line1
        assert len(image2.shape) == 3, err_line1        
    except AssertionError as err:
        print("---\nError")
        print(err)
        
    if vector_of_euler_changes_3D is None:
        print("You must pass in the vector of all possible Euler")
        print("characteristic changes for 3D images, which can be")
        print("computed with euchar.utils.vector_all_euler_changes_in_3D_images.")
        return None

    euler_char_surface = cppsurface.images_3d(image1, image2,
                                              vector_of_euler_changes_3D,
                                              max_intensity1,
                                              max_intensity2)

    return np.array(euler_char_surface)


#=================================================

def bifiltration_density(points, n_neighbors, bins1, bins2,
                         param="alpha"):
    """Compute Euler characteristic surface of bifiltration given by
    combining a given parametrized filtration with the filtration
    derived from a density estimate on points.
    
    Parameters
    -----------
    points
        np.ndarray of shape (N, 2)
    n_neighbors
        int, number of nearest neighbors of a point `p` used to estimate 
        the density at `p`.
    bins1, bins2
        iterable of sorted floats, these are used to bin the 
        parametrization values of simplices in the Alpha filtration. 

    Return
    ------
    euler_char_surface
        np.ndarray of integers

    """
    try:
        err_line1 = "`points` parameter needs to be an np.ndarray\n---"
        assert str(type(points)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    dimension = points.shape[1]
    param_vertices = u.estimate_density(points, n_neighbors)
    
    if dimension == 2:
        if param == "alpha":
            simplices, parametrization = f.alpha_filtration_2d(points)
            euler_char_surf = cppsurface.bifiltration_2d(param_vertices,
                                  parametrization, simplices,
                                  len(bins1), len(bins2),
                                  [min(bins1), max(bins1)],
                                  [min(bins2), max(bins2)])
        else:
            print("Invalid `param` parameter. It can be: 'alpha'.")
            return None
    elif dimension == 3:
        if param == "alpha":
            simplices, parametrization = f.alpha_filtration_3d(points)
            euler_char_surf = cppsurface.bifiltration_3d(param_vertices,
                                  parametrization, simplices,
                                  len(bins1), len(bins2),
                                  [min(bins1), max(bins1)],
                                  [min(bins2), max(bins2)])            
        else:
            print("Invalid `param` parameter. It can be: 'alpha'.")
            return None
    else:
        print("Invalid number of dimensions. `points` must be")
        print("two or three dimensional.")
        return None



    
    
    simplices, param = f.alpha_filtration_2d(points)
    surf = cppsurface.bifiltration_2d(param_vertices,
                                      param,
                                      simplices,
                                      len(bins1), len(bins2),
                                      [min(bins1), max(bins1)],
                                      [min(bins2), max(bins2)])
    return np.array(surf)                


#=================================================

def bifiltration(simplices, parametrization1, parametrization2,
                 bins1, bins2):
    """Compute Euler characteristic surface of bifiltration.
    
    Parameters
    -----------
    simplices
        np.ndarray of shape (N, 3) or (N, 4). Each row is a simplex.
        For example [2 4 -1] is the edge (2,4) and [2 4 9] the 
        triangle (2,4,9). Note that the indices of vertices are 
        sorted.
    parametrization1, parametrization2
        np.ndarray of sorted floats
    bins1, bins2
        np.ndarray of sorted floats, these are used to bin the 
        parametrization values of simplices in the Alpha filtration.
    
    Return
    ------
    euler_char_surface
        np.ndarray of integers

    """
    try:
        err_line1 = "parameters need to be np.ndarrays\n---"
        assert str(type(simplices)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(parametrization1)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(parametrization2)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(bins1)) == "<class 'numpy.ndarray'>", err_line1
        assert str(type(bins2)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    try:
        err_line1 = "simplices must be a two dimensional np.ndarray\n---"
        assert len(simplices.shape) == 2,  err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    dimensions_simplices = u.simplices_to_dimensions(simplices)
    euler_char_surf = cppsurface.bifiltration(dimensions_simplices,
                                              parametrization1,
                                              parametrization2,
                                              bins1, bins2)
    return np.array(euler_char_surf)
