"""


    Implementation of the Euler characateristic curves algorithms.

    Author: Gabriele Beltramo

"""

import numpy as np
import euchar.utils as u
import euchar.filtrations as f
import euchar.cppbinding.curve as cppcurve

#=================================================

def image_2D(image, vector_of_euler_changes_2D=None,
             max_intensity=255):
    """Euler characteristic curve of 2D image.
    
    This uses the vector of all possible Euler characteristic changes 
    produced by a pixel insertion in 2D images. This is recomputed 
    every time this function is called if it is not passes as a 
    parameter.

    Parameters
    ----------
    image
        np.ndarray of integers
    vector_of_euler_changes_2D
        list of integers, precomputed Euler characteristic changes
        produced by a single pixel insertion
    max_intensity
        maximum value of any input of the form of image
    
    Return
    ------
    euler_char_curve
        np.array of integers
    """

    try:
        err_line1 = "`image` parameter needs to be an np.ndarray\n---"
        assert str(type(image)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a two dimensional np.ndarray\n---"
        assert len(image.shape) == 2, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    if vector_of_euler_changes_2D is None:
        vector_of_euler_changes_2D = u.vector_all_euler_changes_in_2D_images()

    euler_char_curve = cppcurve.image_2d(image,
                                         vector_of_euler_changes_2D,
                                         max_intensity)

    return np.array(euler_char_curve)

#=================================================

def image_3D(image, vector_of_euler_changes_3D=None, max_intensity=255):
    """Euler characteristic curve of 3D image.
    
    This uses the vector of all possible Euler characteristic changes 
    produced by a voxel insertion in 3D images. This is must be passed
    as a parameter.
    
    Parameters
    ----------
    image
        np.ndarray of integers
    vector_of_euler_changes_3D
        list of integers, precomputed Euler characteristic changes
        produced by a single voxel insertion
    max_intensity
        maximum value of any input of the form of image
    
    Return
    ------
    euler_char_curve
        np.ndarray of integers

    """

    try:
        err_line1 = "`image` parameter needs to be an np.ndarray\n---"
        assert str(type(image)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a three dimensional np.ndarray\n---"
        assert len(image.shape) == 3, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        
    if vector_of_euler_changes_3D is None:
        print("You must pass in the vector of all possible Euler")
        print("characteristic changes for 3D images, which can be")
        print("computed with euchar.utils.vector_all_euler_changes_in_3D_images.")
        return None

    euler_char_curve = cppcurve.image_3d(image,
                                         vector_of_euler_changes_3D,
                                         max_intensity)

    return np.array(euler_char_curve)

#=================================================

def filtration(points, bins, param="alpha"):
    """Compute the Euler characteristic curve of a finite point set.
    
    Parameters
    ----------
    points
        np.ndarray of shape (N,2) or (N,3)
    bins
        iterable of sorted floats, used to bin the 
        parametrization values of simplices in the Alpha filtration.
    param
        string, 'alpha' for Alpha filtration

    Return
    ------
    euler_char_curve
        np.ndarray of integers
    
    """

    try:
        err_line1 = "`points` parameter needs to be an np.ndarray\n---"
        assert str(type(points)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    dimension = points.shape[1]

    if dimension == 2:
        if param == "alpha":
            simplices, parametrization = f.alpha_filtration_2d(points)
            euler_char_curve = cppcurve.filtration_2d(simplices, parametrization, bins)
        else:
            print("Invalid `param` parameter. It can be: 'alpha'.")
            return None
    elif dimension == 3:
        if param == "alpha":
            simplices, parametrization = f.alpha_filtration_3d(points)
            euler_char_curve = cppcurve.filtration_3d(simplices, parametrization, bins)
        else:
            print("Invalid `param` parameter. It can be: 'alpha'.")
            return None
    else:
        print("Invalid number of dimensions. `points` must be two or")
        print("three dimensional.")
        return None

    return  np.array(euler_char_curve)
