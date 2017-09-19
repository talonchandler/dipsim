import numpy as np

class Fluorophore:
    """A single fluorophore is specified by its 3D position, (unit) absorption
    dipole moment (theta, phi), and (unit) emission dipole moment (theta, phi),
    and ensemble distribution (kappa = None means single molecule. kappa != None
    means Watson distribution).

    """
    def __init__(self, position=np.array([0, 0, 0]),
                 mu_abs=np.array([0, 0]),
                 mu_em=np.array([0, 0]),
                 kappa=None):
        self.position = position
        self.mu_abs = mu_abs        
        self.mu_em = mu_em
        self.mu_ind = 0
        self.kappa = kappa
