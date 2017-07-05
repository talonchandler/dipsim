import numpy as np
import dipsim.util as util

class Detector:
    """A detection path is specified by its optical axis, numerical aperture, and
    index of refraction of the medium.

    """
    def __init__(self, theta_optical_axis, det_type='lens', na=0.8, n=1.33):
        self.theta_optical_axis = theta_optical_axis
        self.det_type = det_type
        self.na = na # Numerical aperture
        self.n = n
        self.f = 10 # Arbitrary
        self.alpha = np.arcsin(self.na/self.n)        
        self.det_rad = self.f*np.tan(self.alpha)

    def calc_collection_efficiency(self, fluorophore):
        if self.det_type=='lens':
            A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
            B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
            theta = util.theta_prime(fluorophore.mu_em[0], fluorophore.mu_em[1], self.theta_optical_axis)
            return 2*(A + B*(np.sin(theta)**2))
        elif self.det_type=='4pi':
            return 1.0
