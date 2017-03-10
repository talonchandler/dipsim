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
        # Populate fisher information matrix
        if self.dist == 'poisson':
            f = self.ev_func
            n = len(x0)
            xx, yy = np.meshgrid(np.arange(n), np.arange(n))
            f_ind = np.stack((xx, yy))

            def calc_fisher_ij(ind, self):
                dx0 = np.eye(1, n, k=ind[0]).flatten()*dx
                dx1 = np.eye(1, n, k=ind[1]).flatten()*dx
                deriv0 = (f(x0+0.5*dx0) - f(x0-0.5*dx0))/np.sum(dx0)
                deriv1 = (f(x0+0.5*dx1) - f(x0-0.5*dx1))/np.sum(dx1)
                return np.sum(deriv0*deriv1/f(x0))
            
            f_inf = np.apply_along_axis(calc_fisher_ij, 0, f_ind, self)
            return np.diag(np.linalg.pinv(f_inf))
