import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, 
    laser, scanned), optical axis, back focal plane source radius, back focal 
    plane polarization, and back focal plane apodization function.

    *** Back focal plane polarization is specified in the x-y plane then 
    rotated. Be careful with oblique illumination. ***

    """
    def __init__(self, illum_type, optical_axis, na=0.8, n=1.33,
                 bfp_pol_dir=None, bfp_apod=None, bfp_n=64):
        
        self.illum_type = illum_type
        
        if np.linalg.norm(optical_axis) - 1.0 > 1e-3:
            print("Warning: optical axis is not a unit vector. Normalizing.")
        self.optical_axis = optical_axis/np.linalg.norm(optical_axis)

        self.na = na
        self.n = n
        self.alpha = np.arcsin(self.na/self.n)
        self.f = 10 # Arbitrary
        self.bfp_rad = self.f*np.tan(self.alpha)

        if np.dot(bfp_pol_dir, np.array([0, 0, 1])) != 0:
            print("Warning: polarization must be specified in x-y plane.")
        elif np.linalg.norm(bfp_pol_dir) - 1 >= 1e-10:
            print("Warning: bfp_pol_dir is not a unit vector. Normalizing.")
        self.bfp_pol_dir = bfp_pol_dir/np.linalg.norm(bfp_pol_dir)

        # Rotate polarization state so that it is perp to optical axis
        self.bfp_pol = np.dot(util.rot_map(self.optical_axis), self.bfp_pol_dir)

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

    def calc_excitation_efficiency(self, fluorophore):
        return np.dot(self.illum_basis, util.mu3to6(fluorophore.mu_abs))
