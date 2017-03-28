from dipsim import fluorophore, illuminator, detector, microscope, stats, util
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
    def __init__(self, microscopes, ):
        self.microscopes = microscopes

    def calc_total_intensities(self, fluorophores):
        I = []
        for m in self.microscopes:
            I.append(m.calc_total_intensity(fluorophores))
        return np.array(I)

    def calc_total_intensities_from_single_fluorophore(self, arguments):
        I = []
        for m in self.microscopes:
            I.append(m.calc_total_intensity_from_single_fluorophore(arguments))
        return np.array(I)
    
    def calc_orientation_std(self, arguments):
        rv = stats.RandomVariable(self.calc_total_intensities_from_single_fluorophore, dist='poisson')
        crlb = rv.crlb(arguments, [1e-2, 1e-2])
        theta = arguments[0]
        theta_std = crlb[0]
        phi_std = crlb[1]
        solid_angle_std = np.sin(theta)*np.sqrt(theta_std)*np.sqrt(phi_std)
        
        return theta_std, phi_std, solid_angle_std

    def plot_orientation_std(self, filename, n=50, my_axs=None, my_caxs=None, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        std_out = np.apply_along_axis(self.calc_orientation_std, 1, directions)
        theta_std = std_out[:,0]
        phi_std = std_out[:,1]
        omega_std = std_out[:,2]
        
        print('Plotting data for microscope: '+filename)
        util.plot_sphere(filename+'_theta.png', directions=directions, data=theta_std, my_ax=my_axs[0], my_cax=my_caxs[0], **kwargs)
        util.plot_sphere(filename+'_phi.png', directions=directions, data=phi_std, my_ax=my_axs[1], my_cax=my_caxs[1], **kwargs)
        util.plot_sphere(filename+'_omega.png', directions=directions, data=omega_std, my_ax=my_axs[2], my_cax=my_caxs[2], **kwargs)

class NFramePolScope(MultiFrameMicroscope):
    """
    An NFramePolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations. 

    An NFramePolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4):
        self.n_frames = n_frames

        # Constant detection path.
        det = detector.Detector(optical_axis=np.array([0, 0, 1]),
                                        na=1.5,
                                        n=1.5)
        # Changing illumination.
        m = []
        for n in range(self.n_frames):
            theta = np.pi*n/self.n_frames
            bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
            ill = illuminator.Illuminator(illum_type='kohler',
                                          optical_axis=np.array([0., 0., 1.]),
                                          f=10,
                                          bfp_rad=3,
                                          bfp_pol=bfp_pol,
                                          bfp_n=4)
            m.append(microscope.Microscope(illuminator=ill, detector=det))
                     
        self.microscopes = m
