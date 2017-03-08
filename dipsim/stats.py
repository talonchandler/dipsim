import numpy as np

class Pdf:
    """
    A probability density function class used to calculate CRLB and other 
    statistics. 
    """
    def __init__(self, func, dist_type='poisson'):
        self.func = func
        self.dist_type = dist_type

    def crlb(self, x0, dx):
        """
        Calculates the Cramer-Rao lower bound at the point x0 using a central 
        difference with width dx.
        
        The length of x0 and dx must match the number of arguments of self.func.
        """
        # Populate fisher information matrix
        if self.dist_type == 'poisson':
            f = self.func
            n = len(x0)
            xx, yy = np.meshgrid(np.arange(n), np.arange(n))
            f_ind = np.stack((xx, yy))

            def calc_fisher_ij(ind, self):
                dx0 = np.eye(1, n, k=ind[0]).flatten()*dx
                dx1 = np.eye(1, n, k=ind[1]).flatten()*dx
                deriv0 = (f(x0+0.5*dx0) - f(x0-0.5*dx0))/np.sum(dx0)
                deriv1 = (f(x0+0.5*dx1) - f(x0-0.5*dx1))/np.sum(dx1)
                return deriv0*deriv1/f(x0)
            
            f_inf = np.apply_along_axis(calc_fisher_ij, 0, f_ind, self)
            return np.diag(np.linalg.pinv(f_inf))
