import numpy as np
import matplotlib.pyplot as plt
import functools

class Microscope:
    """
    A Microscope represents an experiment that collects a single frame of 
    intensity data.  

    A Microscope is specified by its illumination path (an Illuminator object),
    and its detection path (a Detector object).
    """
    def __init__(self, illuminator, detector):
        self.illuminator = illuminator
        self.detector = detector

    def calc_induced_dipoles(self, fluorophores):
        for f in fluorophores:
            f.mu_ind = f.mu_em*np.dot(f.mu_abs, self.illuminator.E_eff)
        
    def calc_total_intensity(self, fluorophores):
        self.calc_induced_dipoles(fluorophores)
        # TODO Green's tensor integrated over area
        # For now sum over entire volume
        I = 0
        for f in fluorophores:
            I += np.linalg.norm(f.mu_ind)**2
        
        return I
        
