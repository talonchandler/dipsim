from dipsim import fluorophore, illuminator, detector, microscope, stats, util, multiframe
import numpy as np
import matplotlib.pyplot as plt
import functools

class MultiFrameMicroscope:
    """A MultiFrameMicroscope represents an experiment that collects multiple 
    frames of intensity data under different conditions (different polarization 
    states or illumination schemes). 

    A MultiFrameMicroscope consists mainly of a list of Microscopes and a
    NoiseModel.
    """
    def __init__(self, ill_thetas, det_thetas, ill_nas, det_nas, ill_types,
                 det_types, colors, n_frames=4, n_pts=1000, max_photons=1000,
                 n_samp=1.33, **kwargs):
        self.n_frames = n_frames
        self.n_pts = n_pts

        m = []
        for i, det_theta in enumerate(det_thetas):
            # Specify detection
            det = detector.Detector(theta_optical_axis=det_thetas[i],
                                    det_type=det_types[i],
                                    na=det_nas[i], n=n_samp)

            # Cycle through polarizations
            for n in range(self.n_frames):
                phi = np.pi*n/self.n_frames
                ill = illuminator.Illuminator(theta_optical_axis=ill_thetas[i],
                                              illum_type=ill_types[i],
                                              na=ill_nas[i], n=n_samp,
                                              phi_pol=phi)

                if type(max_photons) == list:
                    max_photon = max_photons[i]
                else:
                    max_photon = max_photons
                    
                m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photon, color=colors[i]))

        self.microscopes = m
        self.noise_model = stats.NoiseModel(self.calc_total_intensities, **kwargs)
        
    def calc_total_intensities(self, arguments):
        I = []
        for m in self.microscopes:
            I.append(m.calc_intensity(arguments))
        return np.array(I)
    
    def calc_estimation_stats(self):
        sphere_dx = np.arccos(1 - 2/self.n_pts) # avg. half angle betweeen n points on sphere
        self.directions = util.fibonacci_sphere(self.n_pts)

        # Fisher information matrix
        fi_list = []
        for direction in self.directions:
            fi_list.append(self.noise_model.calc_fi(direction, 2*[sphere_dx]))
        self.fi = np.array(fi_list)

        # Solid angle uncertainty
        det = self.fi[:,0,0]*self.fi[:,1,1] - self.fi[:,0,1]*self.fi[:,1,0]
        self.sa_uncert = np.sin(self.directions[:,0])/np.sqrt(det)

        # Coefficient of variation
        self.coeff_of_variation = util.coeff_of_variation(self.sa_uncert)
        
    def scene_string(self):
        asy_string = ''
        for m in self.microscopes:
            asy_string += m.scene_string()
        return asy_string
