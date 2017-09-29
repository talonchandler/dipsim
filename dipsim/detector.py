import numpy as np
import dipsim.util as util

class Detector:
    """A Detector is specified by its optical axis, numerical aperture, and
       index of refraction of the medium.
    """
    def __init__(self, theta_optical_axis=0, det_type='lens', na=0.8,
                 n=1.33, phi_pol=0):
        self.theta_optical_axis = theta_optical_axis
        self.det_type = det_type
        self.na = na # Numerical aperture
        self.n = n
        self.alpha = np.arcsin(self.na/self.n)
        self.phi_pol = phi_pol

    def calc_collection_efficiency(self, fluorophore):
        theta = util.theta_prime(fluorophore.mu_em[0], fluorophore.mu_em[1], self.theta_optical_axis)
        phi = util.phi_prime(fluorophore.mu_em[0], fluorophore.mu_em[1], self.theta_optical_axis)        
        A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
        B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
        C = (7.0/32.0) - (3.0/32.0)*np.cos(self.alpha) - (3.0/32.0)*(np.cos(self.alpha)**2) - (1.0/32.0)*(np.cos(self.alpha)**3)
        if self.det_type=='lens':
            return 2*(A + B*(np.sin(theta)**2))
        elif self.det_type=='polarized':
            return A + B*(np.sin(theta)**2) + C*(np.sin(theta)**2)*np.cos(2*(phi - self.phi_pol))
        elif self.det_type=='4pi':
            return 1.0

