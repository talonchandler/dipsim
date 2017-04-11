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
        theta_std = crlb[0]
        phi_std = crlb[1]
        solid_angle_std = np.sqrt(theta_std)*np.sqrt(phi_std)*np.sin(arguments[0])
        
        return theta_std, phi_std, solid_angle_std

    def plot_orientation_std(self, filename='out.png', n=50, my_ax=None, my_cax=None, **kwargs):
        directions = util.fibonacci_sphere(n)

        std_out = np.apply_along_axis(self.calc_orientation_std, 1, directions, n)
        theta_std = std_out[:,0]
        phi_std = std_out[:,1]
        omega_std = std_out[:,2]

        util.plot_sphere(filename+'_omega.png', directions=directions, data=omega_std, my_ax=my_ax, my_cax=my_cax, **kwargs)

        # Potentially useful for plotting later
        # for i, m in enumerate(self.microscopes):
        #     m.plot_intensities_from_single_fluorophore(str(i)+filename, n, save_file=True)
        # util.plot_sphere(filename+'_theta.png', directions=directions, data=theta_std, my_ax=my_axs[0], my_cax=my_caxs[0], **kwargs)
        # util.plot_sphere(filename+'_phi.png', directions=directions, data=phi_std, my_ax=my_axs[1], my_cax=my_caxs[1], **kwargs)
        
class NFramePolScope(MultiFrameMicroscope):
    """
    An NFramePolScope represents an experiment with the polscope where the
    intensity is collected from N illumination polarizations. 

    An NFramePolScope is specified by the number of frames. 
    """
    def __init__(self, n_frames=4, bfp_n=256, max_photons=100,
                 det_type='lens', det_axis=np.array([0, 0, 1]), **kwargs):
               
        self.n_frames = n_frames

        # Constant detection path.
        det = detector.Detector(optical_axis=det_axis, det_type=det_type, na=1.3, n=1.5)
        
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
                                          bfp_n=bfp_n)
            m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photons))
                     
        MultiFrameMicroscope.__init__(self, m, **kwargs)

    def draw_scene(self, **kwargs):
        pol_dirs = []
        for m in self.microscopes:
            pol_dirs.append(m.illuminator.bfp_pol)
        self.microscopes[0].draw_scene(pol_dirs=pol_dirs, **kwargs)
