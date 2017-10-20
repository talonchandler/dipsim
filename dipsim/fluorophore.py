import numpy as np

class Fluorophore:
    """A fluorophore is specified by its orientation (in theta and phi spherical
    coordinates), it distribution (using a kappa watson distribution), and a
    constant (c) proportional to the fluorohphore's brightness.
    """
    def __init__(self, theta=np.pi/2, phi=0, kappa=None, c=1.0):
        self.theta = theta
        self.phi = phi
        self.kappa = kappa
        self.c = c        
