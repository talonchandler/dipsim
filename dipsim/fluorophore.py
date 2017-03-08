import numpy as np

class Fluorophore:
    """A single fluorophore is specified by its 3D position, (unit) absorption
    dipole moment, and (unit) emission dipole moment.

    """
    def __init__(self, position, mu_abs, mu_em):
        self.position = position
        self.mu_abs = mu_abs        
        self.mu_em = mu_em
        self.mu_ind = 0
