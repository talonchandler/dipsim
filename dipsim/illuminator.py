import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, laser,
    scanned), optical axis, back focal plane source radius, back
    focal plane polarization, and back focal plane apodization function.

    """
    def __init__(self, illum_type, optical_axis, f, bfp_rad, bfp_pol,
                 bfp_apod=None, bfp_n=64):
        
        self.illum_type = illum_type
        
        if np.linalg.norm(optical_axis) != 1.0:
            print("Warning: optical axis is not a unit vector. Normalizing.")
        self.optical_axis = optical_axis/np.linalg.norm(optical_axis)

        if f <= 0:
            print("Warning: f should be positive.")            
        self.f = f

        if bfp_rad <= 0:
            print("Warning: bfp_rad should be positive.")            
        self.bfp_rad = bfp_rad

        if np.dot(bfp_pol, optical_axis) != 0:
            print("Warning: polarization must be orthogonal to optical axis")
        elif np.linalg.norm(bfp_pol) - 1 >= 1e-10:
            print("Warning: bfp_pol is not a unit vector. Normalizing.")
        self.bfp_pol = bfp_pol/np.linalg.norm(bfp_pol)

        # bfp_apod is an apodization function. If no argument is supplied, there
        # is a sharp cutoff at bfp_rad.
        if bfp_apod == None:
            def default_apod(r):
                if r <= self.bfp_rad:
                    return 1.0
                else:
                    return 0
            self.bfp_apod = default_apod
        else:
            self.bfp_apod = bfp_apod

        if bfp_n <= 1:
            print("Warning: bfp_n should be larger than 1. bfp_n should be much\
                   larger than 1 for accurate results. ")
        self.bfp_n = bfp_n
