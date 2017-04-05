import numpy as np
from dipsim import util

class RandomVariable:
    """
    A RandomVariable class used to calculate CRLB and other statistics.

    A RandomVariable is specified by its expected value (which can be a 
    function of other parameters), and its distribution.

    """
    def __init__(self, ev_func, dist='poisson'):
        self.ev_func = ev_func
        self.dist = dist

    def crlb(self, x0, dx, geometry=None):
        """
        Calculates the Cramer-Rao lower bound at the point x0 using a one sided
        difference with width dx. 

        Given the measurement function self.ev_func with noise properties
        self.dist, the CRLB is the minimum variance of an unbiased estimator 
        of the ev_func parameters. 
        
        The length of x0 and dx must match the number of arguments of 
        self.ev_func.

        If two parameters are angles on a sphere, we calculate the derivatives
        differently. Instead of taking derivatives along the coordinate 
        directions, we take derivatives along perpendicular directions on the 
        sphere. 

        The 'geometry' parameter indicates the geometry of each of the 
        input parameters. geometry must have the same length as x0 and dx.

        geometry types:
        'r' = real line parameter (position or intensity)
        't' = theta parameter (the inclination angle)
        'p' = phi parameter (the azimuthal angle)

        Typical inputs if parameter order is [theta, phi, x, y, z] 
        x0 = [0.1, 0.2, 0, 0, 0]
        dx = [1e-2, 1e-2, 1e-2, 1e-2, 1e-2]
        geometry = 'tprrr'
        """
        if self.dist == 'poisson':
            f = self.ev_func
            f0 = f(x0) 
            n = len(x0) # number of params
            m = len(f0) # number of frames (f returns an m x 1 array)

            # Calculate derivative directions
            h = np.zeros((n, n), ) # array of derivative directions
            for i in range(n):
                g = geometry[i]
                if g == 't':
                    ind_t = geometry.index('t')
                    ind_p = geometry.index('p')

                    xyz = util.tp2xyz([x0[ind_t], x0[ind_p]])
                    v1, v2 = util.my_orthonormal_basis(xyz)
                    xyz1 = xyz + v1*dx[ind_t]
                    xyz2 = xyz + v2*dx[ind_p]
                    tp1 = util.xyz2tp(xyz1)
                    tp2 = util.xyz2tp(xyz2)
                    h[ind_t, ind_t] = tp1[ind_t] - x0[ind_t]
                    h[ind_t, ind_p] = tp1[ind_p] - x0[ind_p]
                    h[ind_p, ind_t] = tp2[ind_t] - x0[ind_t]
                    h[ind_p, ind_p] = tp2[ind_p] - x0[ind_p]
                elif g == 'p':
                    break
                elif g == 'r':
                    h[i, :] = np.eye(1, n, k=i).flatten()*dx # ith h vector
                else:
                    print("Warning: geometry must be 't', 'p', or 'r'.")
                    h[i, :] = np.zeros(1, n)
                    
            # Calculate derivative along each direction
            derivs = np.zeros((m, n), )
            for i in range(n):
                derivs[:, i] = (f(x0 + 0.5*h[i, :]) - f(x0 - 0.5*h[i, :]))/dx[i]                
                
            # Calculate fisher information matrix
            f_derivs = np.einsum('ij,ik->ijk', derivs, derivs) # Outer product of each frame
            f_infs = f_derivs/f0[:, np.newaxis, np.newaxis] # Divide each entry by the mean
            f_inf = np.sum(f_infs, axis=0) # Sum over frames

            # Invert and return diagonal
            crlb = np.diag(np.linalg.pinv(f_inf))

            return crlb

