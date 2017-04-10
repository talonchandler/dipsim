import numpy as np

class Detector:
    """A detection path is specified by its optical axis, numerical aperture, and
    index of refraction of the medium.

    """
    def __init__(self, optical_axis, na, n):
        self.optical_axis = optical_axis
        self.na = na # Numerical aperture
        self.n = n

    def calc_collection_efficiency(self, fluorophore):
        alpha = np.arcsin(self.na/self.n)
        A = (1.0/6.0) - (1.0/4.0)*np.cos(alpha) + (1.0/12.0)*(np.cos(alpha)**3)
        B = (1.0/8.0)*np.cos(alpha) - (1.0/8.0)*(np.cos(alpha)**3)
        theta = np.arccos(np.dot(fluorophore.mu_em, self.optical_axis))
        return A + B*(np.sin(theta)**2)
