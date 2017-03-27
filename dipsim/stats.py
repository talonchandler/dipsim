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
            n = len(x0)

            # Evaluate function at all points (2n total evaluations)
            # fpm (n x 2) contains the function at +/- points along each param
            fpm = np.zeros((n, 2))
            for i in range(n):
                v_dx = np.eye(1, n, k=i).flatten()*dx # ith unit vector
                fpm[i, 0] = f(x0 + 0.5*v_dx)
                fpm[i, 1] = f(x0 - 0.5*v_dx)

            # Calculate derivatives
            derivs = (fpm[:,0] - fpm[:,1])/dx

            # Calculate fisher information matrix
            f_inf = np.outer(derivs, derivs)/f(x0)

            # Invert and return diagonal
            crlb =  np.diag(np.linalg.pinv(f_inf))

            return crlb
