import numpy as np
from euchar.utils import vector_all_euler_changes_in_2D_images
from euchar.cppbinding.curve import image_2d, image_3d

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
        list of lists of integers or two dimensional np.array of 
        integers
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
        err_line1 = "`image` parameter needs to be an np.array\n---"
        assert str(type(image)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a two dimensional np.array\n---"
        assert len(image.shape) == 2, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    if vector_of_euler_changes_2D == None:
        vector_of_euler_changes_2D = vector_all_euler_changes_in_2D_images()

    euler_char_curve = image_2d(image, vector_of_euler_changes_2D,
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
        list of lists of list of integers or three dimensional 
        np.array of integers
    vector_of_euler_changes_3D
        list of integers, precomputed Euler characteristic changes
        produced by a single voxel insertion
    max_intensity
        maximum value of any input of the form of image
    
    Return
    ------
    euler_char_curve
        np.array of integers

    """

    try:
        err_line1 = "`image` parameter needs to be an np.array\n---"
        assert str(type(image)) == "<class 'numpy.ndarray'>", err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)

    try:
        err_line1 = "`image` must be a three dimensional np.array\n---"
        assert len(image.shape) == 3, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        
    if vector_of_euler_changes_3d == None:
        print("You must pass in the vector of all possible Euler")
        print("characteristic changes for 3D images, which can be")
        print("computed with euchar.utils.vector_all_euler_changes_in_3D_images.")
        return None

    euler_char_curve = image_3d(image, vector_of_euler_changes_3D,
                                max_intensity)

    return np.array(euler_char_curve)
