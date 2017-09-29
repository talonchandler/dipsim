import numpy as np
import matplotlib.pyplot as plt
import functools
import vispy
from dipsim import util, fluorophore, visuals
from scipy import integrate, special

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

    def calc_intensity(self, direction, kappa=None, epsrel=1e-2):
        fluo = fluorophore.Fluorophore(mu_abs=direction, mu_em=direction)        
        if kappa==None:
            excite = self.illuminator.calc_excitation_efficiency(fluo)
            collect = self.detector.calc_collection_efficiency(fluo)
            return self.max_photons*excite*collect
        else:
            def int_func(theta, phi):
                theta_p = fluo.mu_em[0]
                phi_p = fluo.mu_em[1]
                int_single = self.calc_intensity((theta, phi), kappa=None)
                norm = 1.0/(4*np.pi*special.hyp1f1(0.5, 1.5, kappa))
                weight = np.exp(kappa*(np.dot(util.tp2xyz((theta, phi)), util.tp2xyz((theta_p, phi_p)))**2))
                jacobian = np.sin(theta)
                return jacobian*norm*weight*int_single
            #return integrate.dblquad(int_func, 0, 2*np.pi, lambda x: 0, lambda y: np.pi, epsrel=epsrel, epsabs=0)[0]            
            integral = integrate.nquad(int_func, [[0, np.pi], [0, 2*np.pi]], opts={'epsrel':epsrel, 'epsabs':0, 'limit':1}, full_output=True)
            print(integral)
            return integral[0]

    def plot_intensity(self, filename='', n=50, kappa=None, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        I = np.apply_along_axis(self.calc_intensity, 1, directions, kappa=kappa)
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
            
            pol_string = "arrow(theta, phi_pol, color, false);\n"
            pol_string = pol_string.replace('theta', str(np.rad2deg(ill.theta_optical_axis)))
            pol_string = pol_string.replace('phi_pol', str(np.rad2deg(ill.phi_pol)))
            pol_string = pol_string.replace('color', self.color)        
            asy_string += pol_string
            
        elif ill.illum_type == 'sheet':
            illum_string = "mydot(theta, color);\n"
            illum_string = illum_string.replace('theta', str(ill.theta_optical_axis))
            illum_string = illum_string.replace('color', str(self.color))        
            asy_string += illum_string

            pol_string = "arrow(theta, phi_pol, color, false);\n"
            pol_string = pol_string.replace('theta', str(np.rad2deg(ill.theta_optical_axis)))
            pol_string = pol_string.replace('phi_pol', str(np.rad2deg(ill.phi_pol)))
            pol_string = pol_string.replace('color', self.color)        
            asy_string += pol_string
            
        elif ill.illum_type == 'unpolarized':
            illum_string = "circle(theta, alpha, false, color);\n"
            illum_string = illum_string.replace('theta', str(ill.theta_optical_axis))
            illum_string = illum_string.replace('alpha', str(ill.alpha))
            illum_string = illum_string.replace('color', str(self.color))        
            asy_string += illum_string
            

        # Plot detector
        if det.det_type == 'lens':
            detect_string = "circle(theta, alpha, true, color);\n"
            detect_string = detect_string.replace('theta', str(det.theta_optical_axis))
            detect_string = detect_string.replace('alpha', str(det.alpha+0.01))
            detect_string = detect_string.replace('color', self.color)        
            asy_string += detect_string
        elif det.det_type == 'polarized':
            detect_string = "circle(theta, alpha, true, color);\n"
            detect_string = detect_string.replace('theta', str(det.theta_optical_axis))
            detect_string = detect_string.replace('alpha', str(det.alpha+0.01))
            detect_string = detect_string.replace('color', self.color)        
            asy_string += detect_string
            
            pol_string = "arrow(theta, phi_pol, color, true);\n"
            pol_string = pol_string.replace('theta', str(np.rad2deg(det.theta_optical_axis)))
            pol_string = pol_string.replace('phi_pol', str(np.rad2deg(det.phi_pol)+5))
            pol_string = pol_string.replace('color', self.color)        
            asy_string += pol_string

        sphere_string = """
        // Sphere
        draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));

        // Draw points on sphere
        dotfactor = 7;
        dot(X); 
        dot(Y); 

        circle(0, pi/2, false, (0, 0, 0));
        """

        return asy_string + sphere_string
