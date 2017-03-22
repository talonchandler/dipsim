import numpy as np
import matplotlib.pyplot as plt
import functools
from dipsim import util, fluorophore

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

    def calc_total_intensity_from_single_fluorophore(self, args):
        theta = args[0]
        phi = args[1]

        flu_dir = np.array([np.sin(theta)*np.cos(phi),
                           np.sin(theta)*np.sin(phi),
                           np.cos(theta)])

        flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                      mu_abs=flu_dir,
                                      mu_em=flu_dir)

        I = self.calc_total_intensity([flu])
        return I

    
    def plot_intensities_from_single_fluorophore(self, filename, title, n=50, display='save', random_colors=False, show_edges=False):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_total_intensity_from_single_fluorophore,
                                      1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, title, directions, I, display=display, random_colors=random_colors, show_edges=show_edges)
