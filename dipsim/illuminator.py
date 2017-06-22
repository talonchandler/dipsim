import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, 
    laser, scanned), optical axis, back focal plane source radius, back focal 
    plane polarization, and back focal plane apodization function.

    *** Back focal plane polarization is specified in the x-y plane then 
    rotated. Be careful with oblique illumination. ***

    """
    def __init__(self, illum_type, optical_axis, na=0.8, n=1.33,
                 bfp_pol_dir=None, bfp_apod=None):
        
        self.illum_type = illum_type
        
        if np.linalg.norm(optical_axis) - 1.0 > 1e-3:
            print("Warning: optical axis is not a unit vector. Normalizing.")
        self.optical_axis = optical_axis/np.linalg.norm(optical_axis)

        self.na = na
        self.n = n
        self.alpha = np.arcsin(self.na/self.n)
        self.f = 10 # Arbitrary
        self.bfp_rad = self.f*np.tan(self.alpha)

        if np.dot(bfp_pol_dir, np.array([0, 0, 1])) != 0:
            print("Warning: polarization must be specified in x-y plane.")
        elif np.linalg.norm(bfp_pol_dir) - 1 >= 1e-10:
            print("Warning: bfp_pol_dir is not a unit vector. Normalizing.")
        self.bfp_pol_dir = bfp_pol_dir/np.linalg.norm(bfp_pol_dir)

        # Rotate polarization state so that it is perp to optical axis
        self.bfp_pol = np.dot(util.rot_map(self.optical_axis), self.bfp_pol_dir)

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

    def calc_excitation_efficiency(self, fluorophore):
        A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
        B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
        C = (7.0/32.0) - (3.0/32.0)*np.cos(self.alpha) - (3.0/32.0)*(np.cos(self.alpha)**2) - (1.0/32.0)*(np.cos(self.alpha)**3)        
        D = 4.0/(3.0*(1.0 - np.cos(self.alpha)))
        theta = np.arccos(np.dot(fluorophore.mu_em, self.optical_axis))
        phi = np.arctan2(fluorophore.mu_em[1], fluorophore.mu_em[0])
        phi_pol = np.arctan2(self.bfp_pol_dir[1], self.bfp_pol_dir[0])
        return D*(A + B*(np.sin(theta)**2) + C*(np.sin(theta)**2)*np.cos(2*(phi - phi_pol)))
