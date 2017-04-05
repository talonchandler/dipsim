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
        self.calc_illum_basis()
        
    def calc_illum_basis(self):
        # Generate orthonormal basis with v0 along optical axis
        v0 = self.optical_axis
        v1, v2 = util.orthonormal_basis(v0)

        # Create cartesian sampling of bfp (n x n x 3)
        n = self.bfp_n
        samp = np.linspace(-self.bfp_rad, self.bfp_rad, n)
        xx, yy = np.meshgrid(samp, samp)
        rp = np.einsum('ij,k->ijk', xx, v1) + np.einsum('ij,k->ijk', yy, v2)

        # Find |mu_ind| for each point in bfp            
        def ill_basis_from_bfp_point(rp, ill):
            # Find plane wave normal in front focal plane
            s = ill.optical_axis
            sp = ill.f*s - rp 

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
            apod = ill.bfp_apod(len_rp)

            # Find the rotated GJV
            gjv =  apod*np.dot(R, ill.bfp_pol)

            return util.a3to6(gjv)
        
        ill_basis_rp = np.apply_along_axis(ill_basis_from_bfp_point, 2, rp, self)
        da = (2*self.bfp_rad/n)**2
        a = np.pi*(self.bfp_rad**2)
        self.illum_basis = np.sum(ill_basis_rp, axis=(0, 1))*da/a # Integrate over bfp
