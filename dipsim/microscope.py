import numpy as np
import matplotlib.pyplot as plt
import functools
import vispy
from dipsim import util, visuals, detector, illuminator
from dipsim import fluorophore as flu
from scipy import integrate, special, interpolate
from dipsim.cythonlib.calc_efficiencies import intensity, intensity_from_dist

class Microscope:
    """
    A Microscope represents an experiment that collects a single frame of 
    intensity data.  

    A Microscope is specified by its illumination path (an Illuminator object),
    and its detection path (a Detector object).
    """
    def __init__(self, illuminator=illuminator.Illuminator(),
                 detector=detector.Detector(),
                 max_photons=1, color='(1,0,0)'):
        self.illuminator = illuminator
        self.detector = detector
        self.max_photons = max_photons
        self.color = color
        
        self.precompute_flag = False
        self.fluos_LUT = None
        self.excite_bspl = None
        self.collect_bspl = None
        self.cython_flag = True
    
    def calc_intensity(self, fluorophore=flu.Fluorophore(), epsrel=1e-2):
        # Calculate the intensity measured by this microscope when a given
        # fluorophore distribution is present.
        if np.isinf(fluorophore.kappa):
            if self.cython_flag:
                excite = self.illuminator.calc_excitation_efficiency(fluorophore)
                collect = self.detector.calc_collection_efficiency(fluorophore)
                py_int = self.max_photons*fluorophore.c*excite*collect
                
                c_int = intensity(fluorophore.theta, fluorophore.phi,
                                 self.illuminator.phi_pol,
                                 self.illuminator.theta_optical_axis,                             
                                 self.detector.theta_optical_axis,
                                 self.detector.alpha,
                                 self.max_photons*fluorophore.c)
                return c_int
            else:
                excite = self.illuminator.calc_excitation_efficiency(fluorophore)
                collect = self.detector.calc_collection_efficiency(fluorophore)
                return self.max_photons*fluorophore.c*excite*collect
        else:
            if self.cython_flag:
                c_int = intensity_from_dist(self.illuminator.phi_pol,
                                            self.illuminator.theta_optical_axis,
                                            self.detector.theta_optical_axis,
                                            self.detector.alpha,
                                            self.max_photons*fluorophore.c,
                                            fluorophore.theta,
                                            fluorophore.phi,
                                            fluorophore.kappa)
                return c_int
            else:
                # For fluorophore distributions (this function calls itself!)
                def int_func(theta, phi):
                    int_single = self.calc_intensity(fluorophore=flu.Fluorophore(theta=theta, phi=phi, kappa=np.inf, c=fluorophore.c))
                    norm = 1.0/(4*np.pi*special.hyp1f1(0.5, 1.5, fluorophore.kappa))
                    weight = np.exp(fluorophore.kappa*(np.dot(util.tp2xyz(theta, phi), util.tp2xyz(fluorophore.theta, fluorophore.phi))**2))
                    jacobian = np.sin(theta)
                    return jacobian*weight*int_single
                integral = integrate.nquad(int_func, [[0, np.pi], [0, 2*np.pi]], opts={'epsrel':epsrel, 'epsabs':0, 'limit':1}, full_output=True)
                return integral[0]

    def precompute_efficiencies(self, n_pts=1e3):
        # Generate a LUT for efficiencies
        print("Precomputing efficiencies...")
        directions = util.fibonacci_sphere(n_pts)
        self.fluos = [flu.Fluorophore(*x) for x in directions]
        excite_effs = [self.illuminator.calc_excitation_efficiency(fluo) for fluo in self.fluos]
        collect_effs = [self.detector.calc_collection_efficiency(fluo) for fluo in self.fluos]

        self.excite_bspl = interpolate.bisplrep(directions[:,0], directions[:,1], excite_effs, s=0)
        self.collect_bspl = interpolate.bisplrep(directions[:,0], directions[:,1], collect_effs, s=0)

        self.precompute_flag = True

    def efficiency_from_LUT(self, fluo):
        # Find the closest entry in the lookup table and return
        return interpolate.bisplev(fluo.theta, fluo.phi, self.excite_bspl), interpolate.bisplev(fluo.theta, fluo.phi, self.collect_bspl)
    
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
