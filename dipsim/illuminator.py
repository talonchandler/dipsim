import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, 
    laser, scanned), optical axis, back focal plane source radius, back focal 
    plane polarization, and back focal plane apodization function.

    *** Back focal plane polarization is specified in the x-y plane then 
    rotated. Be careful with oblique illumination. ***

    """
    def __init__(self, illum_type='wide', theta_optical_axis=0, na=0.8, n=1.33, phi_pol=0):
        self.illum_type = illum_type
        self.theta_optical_axis = theta_optical_axis
        self.na = na
        self.n = n
        self.alpha = np.arcsin(self.na/self.n)
        self.f = 10 # Arbitrary
        self.bfp_rad = self.f*np.tan(self.alpha)
        self.phi_pol = phi_pol

    def calc_excitation_efficiency(self, fluorophore):
        theta = util.theta_prime(fluorophore.mu_abs[0], fluorophore.mu_abs[1], self.theta_optical_axis)
        phi = util.phi_prime(fluorophore.mu_abs[0], fluorophore.mu_abs[1], self.theta_optical_axis)
        if self.illum_type == 'wide':
            A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
            B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
            C = (7.0/32.0) - (3.0/32.0)*np.cos(self.alpha) - (3.0/32.0)*(np.cos(self.alpha)**2) - (1.0/32.0)*(np.cos(self.alpha)**3)        
            D = 4.0/(3.0*(1.0 - np.cos(self.alpha)))
            return D*(A + B*(np.sin(theta)**2) + C*(np.sin(theta)**2)*np.cos(2*(phi - self.phi_pol)))
        elif self.illum_type == 'sheet':
            return (np.sin(theta)*np.cos(phi - self.phi_pol))**2
        else:
            print("Warning: invalid illum_type")
