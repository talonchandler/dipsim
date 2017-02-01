import numpy as np
from dipsim import util

class Detector:
    """ A detector class."""
    
    def __init__(self, k, na, n0, n1, optic_axis, img_axis, f_obj, n_pixel, d_pixel):
        self.k = k
        self.na = na # Numerical aperture
        self.n0 = n0 # Glass
        self.n1 = n1 # Medium
        
        print(np.dot(optic_axis, img_axis))
        if np.dot(optic_axis, img_axis) != 0:
            raise ValueError("optic_axis must be orthogonal to img_axis")
        
        self.optic_axis = util.normalize(optic_axis)
        self.img_axis = util.normalize(img_axis)

        self.f_obj = f_obj
        self.n_pixel = n_pixel
        self.d_pixel = d_pixel

    
