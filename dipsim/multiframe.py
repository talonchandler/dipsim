from dipsim import fluorophore, illuminator, detector, microscope, stats, util, multiframe
import numpy as np
import matplotlib.pyplot as plt
import functools

class MultiFrameMicroscope:
    """
    A MultiFrameMicroscope represents an experiment that collects multiple 
    frames of intensity data under different conditions (different polarization 
    states or illumination schemes). 

    A MultiFrameMicroscope consists of a list of Microscopes. 
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
        self.root_det_sin = np.sqrt((self.fi_inv[:,0]*self.fi_inv[:,3] - self.fi_inv[:,1]*self.fi_inv[:,2]))*np.sin(self.directions[:,0])
        
class OneViewPolScope(MultiFrameMicroscope):
    """
    An OneViewPolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations along a single
    illumination and detection arm. 

    An OneViewPolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4, max_photons=1000,
                 illum_type='wide', det_type='lens', illum_det_angle=0,
                 na1=0.8, na2=0.8, n_samp=1.33, **kwargs):
               
        self.n_frames = n_frames

        # Constant detection path.
        ill_axis = illum_det_angle/2
        det_axis = -illum_det_angle/2
        
        det = detector.Detector(theta_optical_axis=det_axis, det_type=det_type,
                                na=na1, n=n_samp)
        
        # Changing illumination.
        m = []
        for n in range(self.n_frames):
            phi = np.pi*n/self.n_frames
            ill = illuminator.Illuminator(illum_type=illum_type,
                                          theta_optical_axis=ill_axis,
                                          na=na2, n=n_samp,
                                          phi_pol=theta)
                                          
            m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons))
                     
        MultiFrameMicroscope.__init__(self, m, **kwargs)
        
    def scene_string(self):
        asy_string = ''
        for m in self.microscopes:
            asy_string += m.scene_string()
        return asy_string

class TwoViewPolScope(MultiFrameMicroscope):
    """
    An MViewPolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations along M 
    """
    def __init__(self, n_frames=4, max_photons=500,
                 illum_type='wide', det_type='lens', illum_det_angle=0,
                 na1=0.8, na2=0.8, n_samp=1.33, **kwargs):
               
        self.n_frames = n_frames

        # Empty list of microscopes
        m = []
        
        # Specify detection and illumination axes.
        axis1 = illum_det_angle/2
        axis2 = -illum_det_angle/2

        # Cycle through detection and illumination orientations
        det_axes = [axis1, axis2]
        det_nas = [na1, na2]
        ill_axes = [axis2, axis1]
        ill_nas = [na2, na1]

        colors = ['(1,0,0)','(0,0,1)']
        
        for det_axis, det_na, ill_axis, ill_na, color in zip(det_axes, det_nas, ill_axes, ill_nas, colors):
            # Specify detection
            det = detector.Detector(theta_optical_axis=det_axis,
                                    det_type=det_type,
                                    na=det_na, n=n_samp)

            # Cycle through polarizations
            for n in range(self.n_frames):
                phi = np.pi*n/self.n_frames
                ill = illuminator.Illuminator(illum_type=illum_type,
                                              theta_optical_axis=ill_axis,
                                              na=ill_na, n=n_samp,
                                              phi_pol=phi)
                
                m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons, color=color))

        MultiFrameMicroscope.__init__(self, m, **kwargs)

    def scene_string(self):
        asy_string = ''
        for m in self.microscopes:
            asy_string += m.scene_string()
        return asy_string
        
