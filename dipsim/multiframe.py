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
    def __init__(self, microscopes, **kwargs):
        self.microscopes = microscopes
        self.noise_model = stats.NoiseModel(self.calc_total_intensities, **kwargs)
        
    def calc_total_intensities(self, arguments):
        I = []
        for m in self.microscopes:
            I.append(m.calc_intensity(arguments))
        return np.array(I)
    
    def calc_orientation_std(self, arguments, n):
        sphere_dx = np.arccos(1 - 2/n) # avg. half angle betweeen n points on sphere
        crlb = self.noise_model.crlb(arguments, [sphere_dx, sphere_dx], geometry='rr')
        return self.noise_model.solid_angle_std(crlb, arguments)

    def plot_orientation_std(self, filename='out.png', n=50, my_ax=None, my_cax=None, **kwargs):
        directions = util.fibonacci_sphere(n)
        omega_std = np.apply_along_axis(self.calc_orientation_std, 1, directions, n)
        util.plot_sphere(filename+'_omega.png', directions=directions, data=omega_std, my_ax=my_ax, my_cax=my_cax, **kwargs)
        
class OneArmPolScope(MultiFrameMicroscope):
    """
    An OneArmPolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations along a single
    illumination and detection arm. 

    An OneArmPolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4, bfp_n=256, max_photons=100,
                 det_type='lens', illum_det_angle=0,
                 f=10, bfp_rad=3, **kwargs):
               
        self.n_frames = n_frames

        # Constant detection path.
        ill_axis = np.array([np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        det_axis = np.array([-np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        
        det = detector.Detector(optical_axis=det_axis, det_type=det_type,
                                na=1.3, n=1.5)
        
        # Changing illumination.
        m = []
        for n in range(self.n_frames):
            theta = np.pi*n/self.n_frames
            bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
            ill = illuminator.Illuminator(illum_type='kohler',
                                          optical_axis=ill_axis,
                                          f=f,
                                          bfp_rad=bfp_rad,
                                          bfp_pol_dir=bfp_pol,
                                          bfp_n=bfp_n)
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
                 f=10, bfp_rad=3, **kwargs):
               
        self.n_frames = n_frames

        # Empty list of microscopes
        m = []
        
        # Specify detection and illumination axes.
        axis1 = np.array([np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])
        axis2 = np.array([-np.sin(illum_det_angle/2), 0, np.cos(illum_det_angle/2)])

        # Cycle through detection and illumination orientations
        det_axes = [axis1, axis2]
        ill_axes = [axis2, axis1]

        for det_axis, ill_axis in zip(det_axes, ill_axes):
            # Specify detection
            det = detector.Detector(optical_axis=det_axis, det_type=det_type,
                                    na=1.3, n=1.5)

            # Cycle through polarizations
            for n in range(self.n_frames):
                theta = np.pi*n/self.n_frames
                bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
                ill = illuminator.Illuminator(illum_type='kohler',
                                              optical_axis=ill_axis,
                                              f=f,
                                              bfp_rad=bfp_rad,
                                              bfp_pol_dir=bfp_pol,
                                              bfp_n=bfp_n)
                m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons))

        MultiFrameMicroscope.__init__(self, m, **kwargs)

    def draw_scene(self, **kwargs):
        pol_dirs = []
        for m in self.microscopes:
            pol_dirs.append(m.illuminator.bfp_pol_dir)
        self.microscopes[0].draw_scene(pol_dirs=pol_dirs, dual_arm=True, **kwargs)
        
