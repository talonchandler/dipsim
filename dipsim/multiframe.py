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

    def calc_solid_angle_min_std(self, args):
        
        def ev_function(self, ev_args):
            theta = ev_args[0]
            phi = ev_args[1]
        
            flu_dir = np.array([np.sin(theta)*np.cos(phi),
                               np.sin(theta)*np.sin(phi),
                               np.cos(theta)])
            
            flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                          mu_abs=flu_dir,
                                          mu_em=flu_dir)
            
            I = self.calc_total_intensities([flu])
            return I
        
        ev_function = functools.partial(ev_function, self)
        rv = stats.RandomVariable(ev_function, dist='poisson')
        crlb = rv.crlb(args, [1e-5, 1e-5])
        theta = args[0]
        theta_std = crlb[0]
        phi_std = crlb[1]
        solid_angle_std = np.sin(theta)*np.sqrt(theta_std)*np.sqrt(phi_std)
        
        return theta_std, phi_std, solid_angle_std

    def plot_solid_angle_min_std(self, filename, n=50, my_axs=None, my_caxs=None, **kwargs):
        directions = util.fibonacci_sphere(n)
        print('Generating data for microscope: '+filename)
        std_out = np.apply_along_axis(self.calc_solid_angle_min_std, 1, directions)
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
