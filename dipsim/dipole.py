import numpy as np

class Dipole:
    """ A dipole class."""
    
    def __init__(self, position, orientation, photon_yield):
        self.position = position
        self.orientation = orientation
        self.photon_yield = photon_yield
