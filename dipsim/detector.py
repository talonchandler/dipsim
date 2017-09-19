import numpy as np
import dipsim.util as util
import dipsim.fluorophore as flu
from scipy import integrate, special

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

    def calc_collection_efficiency(self, fluorophore, epsrel=1e-2):
        if fluorophore.kappa == None:
            theta = util.theta_prime(fluorophore.mu_em[0], fluorophore.mu_em[1], self.theta_optical_axis)
            phi = util.phi_prime(fluorophore.mu_abs[0], fluorophore.mu_abs[1], self.theta_optical_axis)        
            A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
            B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
            C = (7.0/32.0) - (3.0/32.0)*np.cos(self.alpha) - (3.0/32.0)*(np.cos(self.alpha)**2) - (1.0/32.0)*(np.cos(self.alpha)**3)
            if self.det_type=='lens':
                return 2*(A + B*(np.sin(theta)**2))
            elif self.det_type=='polarized':
                return A + B*(np.sin(theta)**2) + C*(np.sin(theta)**2)*np.cos(2*(phi - self.phi_pol))
            elif self.det_type=='4pi':
                return 1.0
        else:
            def int_func(theta, phi):
                theta_p = fluorophore.mu_em[0]
                phi_p = fluorophore.mu_em[1]
                eta_det_single = self.calc_collection_efficiency(flu.Fluorophore(mu_abs=(theta, phi), mu_em=(theta,phi)))
                norm = 1.0/(4*np.pi*special.hyp1f1(0.5, 1.5, fluorophore.kappa))
                weight = np.exp(fluorophore.kappa*(np.dot(util.tp2xyz((theta, phi)), util.tp2xyz((theta_p, phi_p)))**2))
                jacobian = np.sin(theta)
                return jacobian*norm*weight*eta_det_single
            integral = integrate.nquad(int_func, [[0, np.pi], [0, 2*np.pi]], opts={'epsrel':epsrel, 'epsabs':0, 'limit':1}, full_output=True)
            print(integral)
            return integral[0]

