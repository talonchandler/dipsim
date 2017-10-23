import numpy as np
from dipsim import util

class Fluorophore:
    """A fluorophore is specified by its orientation (in theta and phi spherical
    coordinates), it distribution (using a kappa watson distribution), and a
    constant (c) proportional to the fluorohphore's brightness.
    """
    def __init__(self, theta=np.pi/2, phi=0, kappa=np.inf, c=1.0):
        self.theta = theta
        self.phi = phi
        self.kappa = kappa
        self.c = c

    def __sub__(self, x):
        return {'angle_diff': np.rad2deg(util.axis_angle(self.theta, self.phi, x.theta, x.phi)), 
                'kappa_diff': self.kappa - x.kappa,
                'c_diff': self.c - x.c
                }
