import numpy as np

class Detector:
    """A detection path is specified by its optical axis, numerical aperture, and
    index of refraction of the medium.

    """
    def __init__(self, optical_axis, na, n):
        self.optical_axis = optical_axis
        self.na = na # Numerical aperture
        self.n = n
