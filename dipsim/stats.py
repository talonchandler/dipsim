import numpy as np
import scipy.integrate, scipy.misc, scipy.interpolate
from dipsim import util
import functools

class NoiseModel:
    """
    A NoiseModel is specified by an expected value function and the 
    distribution of the measurements. 

    The functions in this class assume that the noise properties are only a
    function of the expected value.
    """
    def __init__(self, ev_func, dist_type='poisson', gauss_mean=0, gauss_std=1,
                 crlb_frame=np.eye(3)):
        self.ev_func = ev_func
        self.dist_type = dist_type

        if dist_type == 'gaussian':
            self.gauss_mean = gauss_mean
            self.gauss_std = gauss_std
        else:
            self.gauss_mean = 0
            self.gauss_std = 0

        self.crlb_frame = crlb_frame
        self.crlb_detector_term_func = self.precompute_detector_term()

    def precompute_detector_term(self):
        if self.dist_type == 'poisson':
            def detector_term_func(mean):
                return 1.0/mean
            
        elif self.dist_type == 'gaussian':
            # Precompute this function for all possible inputs
            out = []
            mean_list = np.array([10.0**x for x in np.arange(-4, 2.1, 0.05)])
            for mean in mean_list:
                def integrand(z, nu, eta, sigma, l_max=250):
                    num = self.pdf(z, nu, eta, sigma, 1, l_max)**2
                    den = self.pdf(z, nu, eta, sigma, 0, l_max)
                    if den != 0:
                        return num/den
                    else:
                        return 1

                min_range = mean-5*(np.sqrt(mean)+self.gauss_std)
                max_range = mean+5*(np.sqrt(mean)+self.gauss_std)                
                integral = scipy.integrate.quad(integrand, min_range, max_range, args=(mean, self.gauss_mean, self.gauss_std))
                out.append(integral[0] - 1)

            # Only perform the integration over the relevant range, otherwise extrapolate
            def whole_detector_term_func(x_in, x_list, y):
                out = []
                for x in x_in:
                    if x < np.max(x_list):
                        f = scipy.interpolate.interp1d(x=x_list, y=y, kind='linear', fill_value='extrapolate')
                        out.append(f(x))
                    else:
                        out.append(1/x)
                return np.array(out)
            
            detector_term_func = functools.partial(whole_detector_term_func, x_list=mean_list, y=out)
            
        return detector_term_func

    def pdf(self, z, nu, eta, sigma, l_min, l_max):
        t1 = np.exp(-nu)/(np.sqrt(2*np.pi)*sigma)
        l = np.arange(l_min, l_max)
        e1 = (l-l_min)*np.log(nu)
        e2 = -np.log(scipy.misc.factorial(l-l_min))
        e3 = -((z-l-eta)**2)/(2*sigma**2)
        out = t1*np.sum(np.exp(e1+e2+e3))
        if np.isfinite(out):
            return out
        else:
            return 0
            
    def plot_detector_term(self, my_ax=None):
        mean_list = np.array([10.0**x for x in np.arange(-6, 4, 0.01)])
        if self.dist_type == 'poisson':
            label = 'Poisson'
        else:
            label = 'P + Gaussian: $\eta$='+str(self.gauss_mean) + ', $\sigma$='+str(self.gauss_std)
        my_ax.plot(mean_list, self.crlb_detector_term_func(mean_list), '-', label=label)
        my_ax.set_xlabel('Mean Number of Electrons')
        my_ax.set_ylabel('Fisher Information Multiplier')        
        my_ax.set_xscale('log')
        my_ax.set_yscale('log')
        my_ax.legend()

    def plot_pdf(self, my_ax=None):
        if self.dist_type == 'poisson':
            label = 'Poisson'
        else:
            label = 'P + Gaussian: $\eta$='+str(self.gauss_mean) + ', $\sigma$='+str(self.gauss_std)
        
        mean = 50
        dx = 0.1
        x = np.arange(mean-4*(np.sqrt(mean)+self.gauss_std), mean+4*(np.sqrt(mean)+self.gauss_std), 0.1)

        pdf_vec = np.vectorize(self.pdf)
        y = pdf_vec(x, mean, self.gauss_mean, self.gauss_std, 0, 100)
        z = pdf_vec(x, mean, self.gauss_mean, self.gauss_std, 1, 100)**2
        
        my_ax.plot(x, y, '-', label=label)
        my_ax.plot(x, z, '-', label=label+'xxx')        
        my_ax.set_xlabel('x')

        my_ax.legend()
        
    def fi_inv(self, x0, dx, geometry=None):
        """
        Calculates the inverse fisher information matrix at the point x0 using 
        a one sided difference with width dx.

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
                tp_prime = util.tp2tp_prime(x0, self.crlb_frame)
                delta = np.eye(1, n, k=i).flatten()*dx # ith h vector
                tp_shifted = util.tp2tp_prime(tp_prime + delta, np.linalg.inv(self.crlb_frame))
                h[i, :] = tp_shifted - x0
            else:
                print("Warning: geometry must be 't', 'p', or 'r'.")
                h[i, :] = np.zeros(1, n)

        # Calculate derivative along each direction
        derivs = np.zeros((m, n), )
        for i in range(n):
            derivs[:, i] = (f(x0 + 0.5*h[i, :]) - f(x0 - 0.5*h[i, :]))/dx[i]                

        # Calculate fisher information matrix
        f_derivs = np.einsum('ij,ik->ijk', derivs, derivs) # Outer product of each frame
        # f_infs = f_derivs/f0[:, np.newaxis, np.newaxis] # Divide each entry by the mean
        f_infs = f_derivs*self.crlb_detector_term_func(f0)[:, np.newaxis, np.newaxis] # Multiply by detector term
        f_inf = np.sum(f_infs, axis=0) # Sum over frames

        # Invert and return
        return np.linalg.pinv(f_inf).flatten()
