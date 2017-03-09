import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, laser,
    scanned), optical axis, back focal plane source radius, back
    focal plane polarization, and back focal plane apodization function.

    """
    def __init__(self, illum_type, optical_axis, f, bfp_rad, bfp_pol,
                 bfp_apod=None, bfp_n=64):
        
        self.illum_type = illum_type
        
        if np.linalg.norm(optical_axis) != 1.0:
            print("Warning: optical axis is not a unit vector. Normalizing.")
        self.optical_axis = optical_axis/np.linalg.norm(optical_axis)

        if f <= 0:
            print("Warning: f should be positive.")            
        self.f = f

        if bfp_rad <= 0:
            print("Warning: bfp_rad should be positive.")            
        self.bfp_rad = bfp_rad

        if np.dot(bfp_pol, optical_axis) != 0:
            print("Warning: polarization must be orthogonal to optical axis")
        elif np.linalg.norm(bfp_pol) - 1 >= 1e-10:
            print("Warning: bfp_pol is not a unit vector. Normalizing.")
        self.bfp_pol = bfp_pol/np.linalg.norm(bfp_pol)

        # bfp_apod is an apodization function. If no argument is supplied, there
        # is a sharp cutoff at bfp_rad.
        if bfp_apod == None:
            def default_apod(r):
                if r <= self.bfp_rad:
                    return 1.0
                else:
                    return 0
            self.bfp_apod = default_apod
        else:
            self.bfp_apod = bfp_apod

        if bfp_n <= 1:
            print("Warning: bfp_n should be larger than 1. bfp_n should be much\
                   larger than 1 for accurate results. ")
        self.bfp_n = bfp_n

        # Calculate E_eff on creation
        self.E_eff = self.calc_E_eff()

    def calc_E_eff(self):
        if self.illum_type == 'kohler':
            # Generate orthonormal basis with v0 along optical axis
            v0 = self.optical_axis
            v1, v2 = util.orthonormal_basis(v0)

            # Create cartesian sampling of bfp (n x n x 3)
            n = self.bfp_n
            samp = np.linspace(-self.bfp_rad, self.bfp_rad, n)
            xx, yy = np.meshgrid(samp, samp)
            rp = np.einsum('ij,k->ijk', xx, v1) + np.einsum('ij,k->ijk', yy, v2)

            # Find E_eff for each point in bfp            
            def E_eff_from_bfp_point(rp, self):
                # Find plane wave normal in ffp
                s = self.optical_axis
                sp = self.f*s - rp 

                # Find rotation matrix
                len_rp = np.linalg.norm(rp)                
                if len_rp == 0:
                    R = np.eye(3)
                else:
                    # Find rotation angle                    
                    theta = np.arccos(np.dot(s, sp/np.linalg.norm(sp)))
                    # Find rotation axis
                    u = np.cross(rp, s)/len_rp 
                    R = util.rot_mat(theta, u) 

                # Find apodization                    
                apod = self.bfp_apod(len_rp)
                # Perform rotation, take abs, and apodize
                return apod*np.abs(np.dot(R, self.bfp_pol)) 

            E_eff_rp = np.apply_along_axis(E_eff_from_bfp_point, 2, rp, self)
            return np.sum(E_eff_rp, axis=(0, 1)) # Sum over bfp
            
        else:
            return np.array([0, 0, 0])
