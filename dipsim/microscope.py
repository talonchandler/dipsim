import numpy as np
import matplotlib.pyplot as plt
from dipsim import util

class Microscope:
    """ A microscope class."""
    
    def __init__(self, exciter, detector):
        self.exciter = exciter
        self.detector = detector
        # Precompute whatever we can

    def calc_bfp(self, n_theta, n_phi, dipole):

        # Create local copies
        d = self.detector 
        s = d.optic_axis
        a = d.img_axis

        # Calculate sample meshgrid (in spherical coordinates)
        phi = np.linspace(0, 2*np.pi, n_phi)
        theta_max = np.arcsin(d.na/d.n0)
        theta = np.linspace(0, theta_max, n_theta)

        phiv, thetav = np.meshgrid(phi, theta)        
        self.phi = phiv
        self.theta = thetav
        
        # Find rotated coordinate axes in terms of cartesian coordinates
        s = self.detector.optic_axis
        x_prime = self.detector.img_axis
        y_prime = np.cross(s, x_prime)

        # Find coordinates of field points r_hat (n_theta x n_phi x 3)
        r_hat = np.einsum('ij,k->ijk', np.cos(thetav), s) + \
                np.einsum('ij,k->ijk', np.cos(phiv)*np.sin(thetav), x_prime) + \
                np.einsum('ij,k->ijk', np.sin(phiv)*np.sin(thetav), y_prime)
        
        # Evaluate the Green's tensor at each field point. "Field of Green's"
        # (n_theta x n_phi x 3 x 3) (thetagrid, phigrid, Green's row, Green's col)

        # Outer product along last dimension
        outer = np.einsum('ijk,ijl->ijkl', r_hat, r_hat) 

        # Create field of unit dyads
        fou = np.tile(np.eye(3), (n_theta, n_phi, 1, 1))

        # Calculate phase term
        rp_vec = dipole.position
        rp_dot_rhat = np.einsum('i,jki->jk', rp_vec, r_hat)
        phase = np.exp(1j*d.k*d.n0*(d.f_obj - rp_dot_rhat))/(4*np.pi*d.f_obj)
        phase_tiled = np.tile(phase, (3, 3, 1, 1)).transpose([2,3,0,1])

        # Create field of Green's
        fog = phase_tiled*(fou - outer)
        
        # Calculate apodization factor
        apod = np.sqrt(d.n1/(d.n0*np.einsum('ijk,k->ij', r_hat, s)))
        apod3 = np.tile(apod, (3, 1, 1)).transpose([1, 2, 0])

        # Calculate rotation matrix
        
        # Multiply Green's matrix by the dipole
        E_bfp = apod3*np.einsum('ijkl,l->ijk', fog, dipole.orientation)
        I_bfp = np.einsum('ijk,ijk->ij', E_bfp, E_bfp.conjugate())

        self.E_bfp = E_bfp
        self.I_bfp = I_bfp.real
        
    def calc_img(self):

        # Sample on cartesian grid
        
        E_img = np.fft.fftshift(np.fft.fftn(self.E_bfp, axes=(0, 1)), axes=(0,1))
        I_img = np.einsum('ijk,ijk->ij', E_img, E_img.conjugate())

        self.E_img = I_img
        self.I_img = I_img.real


    def plot_detector(self, filename):
        import matplotlib.colors as colors
        from matplotlib.colors import BoundaryNorm
        from matplotlib.ticker import MaxNLocator

        f, (ax0, ax1) = plt.subplots(1, 2, figsize=(5,5),
                                     subplot_kw=dict(projection='polar'))
        ax0 = plt.subplot(121, projection='polar')
        ax1 = plt.subplot(122)
        
        ax0.set_xlabel(str(np.round(self.detector.img_axis, 2)))
        ax0.set_ylabel(str(np.round(np.cross(self.detector.optic_axis,
                                             self.detector.img_axis),2)))
        ax0.set_title('Back Focal Plane')
        ax0.grid(False)
        ax0.get_xaxis().set_ticks([])
        ax0.get_yaxis().set_ticks([])

        levels = MaxNLocator(nbins=256).tick_values(self.I_bfp.min(), self.I_bfp.max())        
        cmap = plt.get_cmap('gray')
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
        
        im = ax0.pcolor(self.phi, self.theta, self.I_bfp,
                        norm=norm, cmap=cmap)

        ax1.set_xlabel(str(np.round(self.detector.img_axis, 2)))
        ax1.set_ylabel(str(np.round(np.cross(self.detector.optic_axis,
                                             self.detector.img_axis))))
        ax1.set_title('Detector Plane', y=1.05)
        ax1.grid(False)
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        
        ax1.imshow(np.zeros((100, 100)),
                   origin='lower',
                   interpolation='none')
        
        print("Saving:", filename)
        f.savefig(filename)

