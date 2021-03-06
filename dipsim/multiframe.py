from dipsim import illuminator, detector, microscope, stats, util
import dipsim.fluorophore as flu
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
    def __init__(self, ill_thetas=[0], det_thetas=[0],
                 ill_nas=[0.8], det_nas=[0.8],
                 ill_types=['wide'], det_types=['lens'],
                 colors=['(1,0,0)'], n_frames=4, n_det_frames=1,
                 frame_offsets=None, max_photons=1, n_samp=1.33,
                 **kwargs):
        
        self.n_frames = n_frames
        self.frame_offsets = frame_offsets
        self.n_det_frames = n_det_frames

        m = []
        for i, det_theta in enumerate(det_thetas):
            # Cycle through detection polarizations
            for k in range(self.n_det_frames):
                phi_det = np.pi*k/self.n_det_frames
                det = detector.Detector(theta_optical_axis=det_thetas[i],
                                        det_type=det_types[i],
                                        na=det_nas[i], n=n_samp,
                                        phi_pol=phi_det)

                # Cycle through illumination polarizations
                for n in range(self.n_frames):
                    phi_ill = np.pi*n/self.n_frames
                    if self.frame_offsets != None:
                        if self.frame_offsets[i]:
                            phi_ill += np.pi/(2*self.n_frames)
                            
                    ill = illuminator.Illuminator(theta_optical_axis=ill_thetas[i],
                                                  illum_type=ill_types[i],
                                                  na=ill_nas[i], n=n_samp,
                                                  phi_pol=phi_ill)

                    if type(max_photons) == list:
                        max_photon = max_photons[i]
                    else:
                        max_photon = max_photons

                    m.append(microscope.Microscope(illuminator=ill, detector=det, max_photons=max_photon, color=colors[i]))

        self.microscopes = m
        self.noise_model = stats.NoiseModel(self.calc_intensities, **kwargs)
        
    def calc_intensities(self, fluorophore=flu.Fluorophore(), epsrel=1e-2):
        return [x.calc_intensity(fluorophore=fluorophore, epsrel=epsrel) for x in self.microscopes]
    
    def calc_estimation_stats(self, n_pts):
        sphere_dx = np.arccos(1 - 2/n_pts) # avg. half angle betweeen n points on sphere
        directions = util.fibonacci_sphere(n_pts)

        # Fisher information matrix
        fi_list = []
        for direction in directions:
            fi_list.append(self.noise_model.calc_fi(direction, 2*[sphere_dx]))
        self.fi = np.array(fi_list)

        # Solid angle uncertainty
        det = self.fi[:,0,0]*self.fi[:,1,1] - self.fi[:,0,1]*self.fi[:,1,0]
        self.sa_uncert = np.sin(self.directions[:,0])/np.sqrt(det)

        # Coefficient of variation
        self.coeff_of_variation = util.coeff_of_variation(self.sa_uncert)

    def precompute_efficiencies(self, n_pts=1e3):
        for m in self.microscopes:
            m.precompute_efficiencies(n_pts)

    def scene_string(self):
        asy_string = ''
        for m in self.microscopes:
            asy_string += m.scene_string()
        return asy_string

    def ellipse_string(self, n_pts=20):
        pts = util.fibonacci_sphere(n_pts)
        asy_string = ''
        for pt in pts:
            fi = self.noise_model.calc_fi(pt, 2*[self.sphere_dx])
            fi_inv = np.linalg.pinv(fi)
            if fi_inv[0,0] < 0:
                fi_inv[0,0] = 0
            sigx = np.sqrt(fi_inv[0,0])
            sigxy = np.sin(pt[0])*fi_inv[1,0]
            sigy = np.sin(pt[0])*np.sqrt(fi_inv[1,1])
            a = np.sqrt(0.5*(sigx**2+sigy**2) + np.sqrt(0.25*((sigx**2-sigy**2)**2) + sigxy**2))
            b = np.sqrt(0.5*(sigx**2+sigy**2) - np.sqrt(0.25*((sigx**2-sigy**2)**2) + sigxy**2))
            theta = 0.5*np.arctan2(2*sigxy,sigx**2-sigy**2)
            
            ellipse_string = 'ellipse(Theta, Phi, a1, b1, theta, false, (0,0,0));\n'
            ellipse_string  = ellipse_string.replace('Theta', str(pt[0]))
            ellipse_string  = ellipse_string.replace('Phi', str(pt[1]))
            ellipse_string  = ellipse_string.replace('a1', str(a))
            ellipse_string  = ellipse_string.replace('b1', str(b))
            ellipse_string  = ellipse_string.replace('theta', str(theta))
            ellipse_string  = ellipse_string.replace('grey', str(theta))

            if np.dot(util.tp2xyz(pt), 3*[1.0/np.sqrt(3)]) > 0:
                asy_string += ellipse_string
            
        return asy_string
