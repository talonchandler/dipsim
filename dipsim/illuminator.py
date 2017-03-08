import numpy as np
import dipsim.util as util

class Illuminator:
    """An illumination path is specified by its illumination type (Kohler, laser,
    scanned), optical axis, back focal plane source radius, and back
    focal plane polarizer.

    """
    def __init__(self, illum_type, optical_axis, f, bfp_rad, bfp_pol, bfp_n=64):
        self.illum_type = illum_type
        
        if np.linalg.norm(optical_axis) != 1.0:
            print("Warning: optical axis is not a unit vector. Normalizing.")
        self.optical_axis = optical_axis/np.linalg.norm(optical_axis)
        
        self.f = f
        self.bfp_rad = bfp_rad
        self.bfp_n = bfp_n

        if np.dot(bfp_pol, optical_axis) != 0:
            print("Warning: polarization must be orthogonal to optical axis")
        elif np.linalg.norm(bfp_pol) - 1 >= 1e-10:
            print("Warning: bfp_pol is not a unit vector. Normalizing.")
        self.bfp_pol = bfp_pol/np.linalg.norm(bfp_pol)

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
                s = self.optical_axis
                sp = self.f*s - rp # Find new direction of prop
                theta = np.arccos(np.dot(s, sp/np.linalg.norm(sp))) # Find rotation angle
                u = np.cross(rp, s)/np.linalg.norm(rp) # Find rotation axis
                R = util.rot_mat(theta, u) # Find rotation matrix
                return np.abs(np.dot(R, self.bfp_pol)) # Perform rotation and take abs
                            
            E_eff_rp = np.apply_along_axis(E_eff_from_bfp_point, 2, rp, self)
            E_eff = np.sum(E_eff_rp, axis=(0, 1))

            print(E_eff/np.linalg.norm(E_eff))
            return E_eff/np.linalg.norm(E_eff)
        else:
            return np.array([0, 0, 0])
