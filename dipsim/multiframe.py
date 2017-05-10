from dipsim import fluorophore, illuminator, detector, microscope, stats, util, multiframe
import numpy as np
import matplotlib.pyplot as plt
import functools

class MultiFrameMicroscope:
    """
    A MultiFrameMicroscope represents an experiment that collects multiple 
    frames of intensity data under different conditions (different polarization 
    states or illumination schemes). 

    A MultiFrameMicroscope is specified by a list of Microscopes. 
    """
    def __init__(self, microscopes, n_pts=1000, **kwargs):
        self.microscopes = microscopes
        self.noise_model = stats.NoiseModel(self.calc_total_intensities, **kwargs)
        self.calc_estimation_stats(n_pts)
        
    def calc_total_intensities(self, arguments):
        I = []
        for m in self.microscopes:
            I.append(m.calc_intensity(arguments))
        return np.array(I)
    
    def calc_estimation_stats(self, n):
        def calc_single_fi_inv(arguments, self, n):
            sphere_dx = np.arccos(1 - 2/n) # avg. half angle betweeen n points on sphere
            return self.noise_model.fi_inv(arguments, [sphere_dx, sphere_dx], geometry='rr')

        self.directions = util.fibonacci_sphere(n)
        self.fi_inv = np.apply_along_axis(calc_single_fi_inv, 1, self.directions, self, n)
        self.sa_uncert = np.sqrt(self.fi_inv[:,0])*np.sqrt(self.fi_inv[:,3])*np.sin(self.directions[:,0])
        
class OneArmPolScope(MultiFrameMicroscope):
    """
    An OneArmPolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations along a single
    illumination and detection arm. 

    An OneArmPolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4, bfp_n=256, max_photons=100,
                 det_type='lens', illum_det_angle=0,
                 na_ill=0.8, na_det=0.8, n_samp=1.33, **kwargs):
               
        self.n_frames = n_frames

        # Constant detection path.
        ill_axis = np.array([np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        det_axis = np.array([-np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        
        det = detector.Detector(optical_axis=det_axis, det_type=det_type,
                                na=na_det, n=n_samp)
        
        # Changing illumination.
        m = []
        for n in range(self.n_frames):
            theta = np.pi*n/self.n_frames
            bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
            ill = illuminator.Illuminator(illum_type='kohler',
                                          optical_axis=ill_axis,
                                          na=na_ill, n=n_samp,
                                          bfp_pol_dir=bfp_pol, bfp_n=bfp_n)
                                          
            m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons))
                     
        MultiFrameMicroscope.__init__(self, m, **kwargs)

    def draw_scene(self, **kwargs):
        pol_dirs = []
        for m in self.microscopes:
            pol_dirs.append(m.illuminator.bfp_pol_dir)
        self.microscopes[0].draw_scene(pol_dirs=pol_dirs, **kwargs)

class TwoArmPolScope(MultiFrameMicroscope):
    """
    A TwoArmPolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations along an
    illumination and detection arm then the illumination and detection arms 
    are swapped.

    A TwoArmPolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4, bfp_n=256, max_photons=100,
                 det_type='lens', illum_det_angle=0,
                 na1=0.8, na2=0.8, n_samp=1.33, **kwargs):
               
        self.n_frames = n_frames

        # Empty list of microscopes
        m = []
        
        # Specify detection and illumination axes.
        axis1 = np.array([np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        axis2 = np.array([-np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])

        # Cycle through detection and illumination orientations
        det_axes = [axis1, axis2]
        det_nas = [na1, na2]
        ill_axes = [axis2, axis1]
        ill_nas = [na2, na1]
        
        for det_axis, det_na, ill_axis, ill_na in zip(det_axes, det_nas, ill_axes, ill_nas):
            # Specify detection
            det = detector.Detector(optical_axis=det_axis, det_type=det_type,
                                    na=det_na, n=n_samp)

            # Cycle through polarizations
            for n in range(self.n_frames):
                theta = np.pi*n/self.n_frames
                bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
                ill = illuminator.Illuminator(illum_type='kohler',
                                              optical_axis=ill_axis,
                                              na=ill_na, n=n_samp,
                                              bfp_pol_dir=bfp_pol,
                                              bfp_n=bfp_n)
                
                m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons))

        MultiFrameMicroscope.__init__(self, m, **kwargs)

    def draw_scene(self, **kwargs):
        pol_dirs = []
        for m in self.microscopes:
            pol_dirs.append(m.illuminator.bfp_pol_dir)
        self.microscopes[0].draw_scene(pol_dirs=pol_dirs, dual_arm=True, **kwargs)
        
