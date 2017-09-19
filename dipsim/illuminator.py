import numpy as np
import dipsim.util as util
import dipsim.fluorophore as flu
from scipy import integrate, special

class Illuminator:
    """An Illuminator is specified by its illumination type (wide, sheet),
    optical axis, numerical aperture, index of refraction of the sample, and
    polarization.
    """
    def __init__(self, illum_type='wide', theta_optical_axis=0, na=0.8,
                 n=1.33, phi_pol=0):
        self.illum_type = illum_type
        self.theta_optical_axis = theta_optical_axis
        self.na = na
        self.n = n
        self.alpha = np.arcsin(self.na/self.n)
        self.phi_pol = phi_pol

    def calc_excitation_efficiency(self, fluorophore, epsrel=1e-2):
        if fluorophore.kappa == None:
            theta = util.theta_prime(fluorophore.mu_abs[0], fluorophore.mu_abs[1], self.theta_optical_axis)
            phi = util.phi_prime(fluorophore.mu_abs[0], fluorophore.mu_abs[1], self.theta_optical_axis)
            A = (1.0/4.0) - (3.0/8.0)*np.cos(self.alpha) + (1.0/8.0)*(np.cos(self.alpha)**3)
            B = (3.0/16.0)*np.cos(self.alpha) - (3.0/16.0)*(np.cos(self.alpha)**3)
            C = (7.0/32.0) - (3.0/32.0)*np.cos(self.alpha) - (3.0/32.0)*(np.cos(self.alpha)**2) - (1.0/32.0)*(np.cos(self.alpha)**3)        
            D = 4.0/(3.0*(1.0 - np.cos(self.alpha)))
            if self.illum_type == 'wide':
                return D*(A + B*(np.sin(theta)**2) + C*(np.sin(theta)**2)*np.cos(2*(phi - self.phi_pol)))
            elif self.illum_type == 'sheet':
                return (np.sin(theta)*np.cos(phi - self.phi_pol))**2
            elif self.illum_type == 'unpolarized':
                return 2*D*(A + B*(np.sin(theta)**2))
            else:
                print("Warning: invalid illum_type")
        else:
            def int_func(theta, phi):
                theta_p = fluorophore.mu_abs[0]
                phi_p = fluorophore.mu_abs[1]
                eta_abs_single = self.calc_excitation_efficiency(flu.Fluorophore(mu_abs=(theta, phi), mu_em=(theta,phi)))
                norm = 1.0/(4*np.pi*special.hyp1f1(0.5, 1.5, fluorophore.kappa))
                weight = np.exp(fluorophore.kappa*(np.dot(util.tp2xyz((theta, phi)), util.tp2xyz((theta_p, phi_p)))**2))
                jacobian = np.sin(theta)
                return jacobian*norm*weight*eta_abs_single
            return integrate.dblquad(int_func, 0, 2*np.pi, lambda x: 0, lambda y: np.pi, epsrel=epsrel, epsabs=0)[0]
