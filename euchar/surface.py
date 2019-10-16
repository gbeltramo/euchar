import numpy as np
import euchar.utils as u
import euchar.filtrations as f
import euchar.curve as c
from euchar.cppbinding.utils import char_binary_image_2d
import euchar.cppbinding.surface as cppsurface

#=================================================

def images_2D(image1, image2, vector_of_euler_changes_2D=None,
              max_intensity1=255, max_intensity2=255):
    """
    Euler characteristic curve of a pair of 2D images.
    
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
        return None

    try:
        err_line1 = "`image` must be a two dimensional np.array\n---"
        assert len(image1.shape) == 2, err_line1
        assert len(image2.shape) == 2, err_line1
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    if vector_of_euler_changes_2D is None:
        vector_of_euler_changes_2D = u.vector_all_euler_changes_in_2D_images()

    euler_char_surface = cppsurface.images_2d(image1, image2,
                                            vector_of_euler_changes_2D,
                                            max_intensity1,
                                            max_intensity2)

    return np.array(euler_char_surface)

#=================================================

def image_2D_and_function(image, func, max_intensity=255,
                          iterations=1, kwargs=None):
    """
    Compute the Euler characteristic surface of a give image with the
    bi-filtration given by sublevel sets of pixel values and 
    iterative application of a given function on these sublevel sets.
    
    Parameters
    ----------
    image
        np.ndarray of integers
    iterations
        number of times the function is applied to the image.
    func
        function taking a sublevel set of `image` and returning 
        a binary image of the same size. 
    max_intensity
        maximum value of elements in image
    kwargs
        possible arguments of function `func` passed as a parameter.

    Returns
    -------
    euler_char_surface
        np.ndarray of integers
    
    """
    
    surface = np.zeros(shape=(iterations+1, max_intensity+1), dtype=int)

    surface[0] = c.image_2D(image, max_intensity=max_intensity)
    
    if kwargs is None:
        kwargs = dict()
        
    for i in range(1, iterations+1):
        for j in range(max_intensity+1):
            thresh = (image <= j).astype(np.uint8)
            for _ in range(i):
                thresh = func(thresh, **kwargs)
            euler_char = char_binary_image_2d(thresh.astype(bool))
            surface[i][j] = euler_char
            
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
        maximum values of elements in image1 and image2
    
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
        return None

    try:
        err_line1 = "`image` must be a three dimensional np.array\n---"
        assert len(image1.shape) == 3, err_line1
        assert len(image2.shape) == 3, err_line1        
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    try:
        err_line1 = "You must pass in the vector of all possible Euler"
        err_line2 = "characteristic changes for 3D images, which can be"
        err_line3 = "computed with euchar.utils.vector_all_euler_changes_in_3D_images.\n---"
        assert vector_of_euler_changes_3D is not None, err_line1 + err_line2 + err_line3
    except AssertionError as err:
        print("---\nError")
        print(err)
        return None

    euler_char_surface = cppsurface.images_3d(image1, image2,
                                              vector_of_euler_changes_3D,
                                              max_intensity1,
                                              max_intensity2)

    return np.array(euler_char_surface)

#=================================================

def bifiltration(simplices, parametrization1, parametrization2,
                 bins1, bins2):
    """
    Compute Euler characteristic surface of bifiltration.
    
    Parameters
    -----------
    simplices
        np.ndarray of shape (N, 3) or (N, 4). Each row is a simplex.
        For example [2 4 -1] is the edge (2,4) and [2 4 9] the 
        triangle (2,4,9). 
    parametrization1, parametrization2
        np.ndarray of floats
    bins1, bins2
        np.ndarrays of sorted floats, used to bin the parametrization
        values of simplices 
    
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
