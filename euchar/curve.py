from euchar.cppbinding import utils, curve

#=================================================

def image_2d(image, vector_of_euler_changes_2d=None, max_intensity=256):
    """Euler characteristic curve of image.

    Parameters
    ----------
    image
        list of lists of integers or two dimensional np.array
    vector_of_euler_changes_2d
        list of integers, precomputed Euler characteristic changes
        produced by a single pixel insertion
    max_intensity
        maximum value of any input of the form of image
    
    Return
    ------
    euler_char_curve
        list of integers
    """

    if vector_of_euler_changes_2d == None:
        vector_of_euler_changes_2d = utils.vector_of_euler_changes_2d()

    euler_char_curve = curve.image_2d(image, max_intensity, vector_of_euler_changes_2d)

    return euler_char_curve

#=================================================

def complex_2d():
    pass
