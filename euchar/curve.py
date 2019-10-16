import numpy as np
import euchar.utils as u
import euchar.filtrations as f
import euchar.cppbinding.curve as cppcurve

#=================================================

def image_2D(image, vector_of_euler_changes_2D=None,
             max_intensity=255):
    """
    Euler characteristic curve of 2D image.
    
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
        maximum value of elements in image
    
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
    """
    Euler characteristic curve of 3D image.
    
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
        maximum value of elements in image
    
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
        return None

    try:
        err_line1 = "`image` must be a three dimensional np.ndarray\n---"
        assert len(image.shape) == 3, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    try:
        err_line1 = "You must pass in the vector of all possible Euler\n"
        err_line2 = "characteristic changes for 3D images, which can be\n"
        err_line3 = "computed with euchar.utils.vector_all_euler_changes_in_3D_images.\n---"
        assert vector_of_euler_changes_3D is not None, err_line1 + err_line2 + err_line3
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    euler_char_curve = cppcurve.image_3d(image,
                                         vector_of_euler_changes_3D,
                                         max_intensity)

    return np.array(euler_char_curve)

#=================================================

def filtration(simplices, parametrization, bins):
    """
    Compute the Euler characteristic curve of a filtration.
    
    Parameters
    ----------
    simplices
        np.ndarray of shape (N, 3) or (N, 4). Each row is a simplex.
        For example [2 4 -1] is the edge (2,4) and [2 4 9] the 
        triangle (2,4,9). 
    parametrization
        np.ndarray of floats
    bins
        np.ndarray of sorted floats, used to bin the parametrization
        values of simplices 
    
    Return
    ------
    euler_char_curve
        np.ndarray of integers

    """

    try:
        err_line1 = "`simplices` and `parametrization` must have same\n"
        err_line2 = "length."
        assert len(simplices) == len(parametrization), err_line1+err_line2
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    dim_simplices = u.simplices_to_dimensions(simplices) 
    euler_char_curve = cppcurve.filtration(dim_simplices, parametrization, bins)
        
    return  np.array(euler_char_curve)
