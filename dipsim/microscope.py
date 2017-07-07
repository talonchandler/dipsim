import numpy as np
import matplotlib.pyplot as plt
import functools
import vispy
from dipsim import util, fluorophore, visuals

class Microscope:
    """
    A Microscope represents an experiment that collects a single frame of 
    intensity data.  

    A Microscope is specified by its illumination path (an Illuminator object),
    and its detection path (a Detector object).
    """
    def __init__(self, illuminator, detector, max_photons, color='(1,0,0)'):
        self.illuminator = illuminator
        self.detector = detector
        self.max_photons = max_photons
        self.color = color

    def calc_intensity(self, direction):
        return self.max_photons*self.calc_sensitivity(direction)
    
    def calc_sensitivity(self, direction):
        flu = fluorophore.Fluorophore(mu_abs=direction, mu_em=direction)
        excite = self.illuminator.calc_excitation_efficiency(flu)
        collect = self.detector.calc_collection_efficiency(flu)        
        return excite*collect

    def plot_sensitivity(self, filename='', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_sensitivity, 1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def calc_excitation_efficiency(self, direction):
        flu = fluorophore.Fluorophore(mu_abs=direction, mu_em=direction)
        I = self.illuminator.calc_excitation_efficiency(flu)
        return I
    
    def plot_excitation_efficiency(self, filename='', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_excitation_efficiency, 1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def calc_collection_efficiency(self, direction):
        flu = fluorophore.Fluorophore(mu_abs=direction, mu_em=direction)
        I = self.detector.calc_collection_efficiency(flu)
        return I
    
    def plot_collection_efficiency(self, filename='', n=50, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_collection_efficiency, 1, directions)
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename, directions=directions, data=I, **kwargs)

    def scene_string(self):
        ill = self.illuminator
        det = self.detector
        asy_string = ''
        if ill.illum_type == 'wide':            
            illum_string = "circle(theta, alpha, false, color);\n"
            illum_string = illum_string.replace('theta', str(ill.theta_optical_axis))
            illum_string = illum_string.replace('alpha', str(ill.alpha))
            illum_string = illum_string.replace('color', str(self.color))        
            asy_string += illum_string
        elif ill.illum_type == 'sheet':
            illum_string = "mydot(theta, color);\n"
            illum_string = illum_string.replace('theta', str(ill.theta_optical_axis))
            illum_string = illum_string.replace('color', str(self.color))        
            asy_string += illum_string

        # Plot detector
        detect_string = "circle(theta, alpha, true, color);\n"
        detect_string = detect_string.replace('theta', str(det.theta_optical_axis))
        detect_string = detect_string.replace('alpha', str(det.alpha+0.01))
        detect_string = detect_string.replace('color', self.color)        
        asy_string += detect_string

        # Plot polarization
        pol_string = "arrow(theta, phi_pol, color);\n"
        pol_string = pol_string.replace('theta', str(np.rad2deg(ill.theta_optical_axis)))
        pol_string = pol_string.replace('phi_pol', str(np.rad2deg(ill.phi_pol)))
        pol_string = pol_string.replace('color', self.color)        
        asy_string += pol_string

        return asy_string
