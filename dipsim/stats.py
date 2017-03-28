import numpy as np

class RandomVariable:
    """
    A RandomVariable class used to calculate CRLB and other statistics.

    A RandomVariable is specified by its expected value (which can be a 
    function of other parameters), and its distribution.

    """
    def __init__(self, ev_func, dist='poisson'):
        self.ev_func = ev_func
        self.dist = dist

    def crlb(self, x0, dx):
        """
        Calculates the Cramer-Rao lower bound at the point x0 using a central
        difference with width dx. 

        Given the measurement function self.ev_func with noise properties
        self.dist, the CRLB is the minimum variance of an unbiased estimator 
        of the ev_func parameters. 
        
        The length of x0 and dx must match the number of arguments of 
        self.ev_func.
        """
        if self.dist == 'poisson':
            f = self.ev_func
            f0 = f(x0) 
            n = len(x0) # number of params
            m = len(f0) # number of frames (f returns an m x 1 array)

            # Calculate the derivative along each direction
            derivs = np.zeros((m, n), )
            for i in range(n):
                h = np.eye(1, n, k=i).flatten()*dx # ith h vector
                derivs[:, i] = (f(x0 + h) - f0)/h[i]

            # Calculate fisher information matrix
            f_derivs = np.einsum('ij,ik->ijk', derivs, derivs) # Outer product of each frame
            f_infs = f_derivs/f0[:, np.newaxis, np.newaxis] # Divide each entry by the mean
            f_inf = np.sum(f_infs, axis=0) # Sum over frames

            # Invert and return diagonal
            crlb = np.diag(np.linalg.pinv(f_inf))

            return crlb
