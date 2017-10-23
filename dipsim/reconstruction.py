import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from dipsim import util, multiframe, fluorophore

class Reconstruction():
    """A reconstruction is specified by a MuliFrameMicroscope, data, and a 
    reconstruction options. 
    """
    def __init__(self, multiframe=multiframe.MultiFrameMicroscope(),
                 data=None, recon_type=None, recon_dx=1e-6, recon_eps=1e-3,
                 recon_max_iter=1e2, recon_start=None):
        
        self.multiframe = multiframe
        self.data = data
        if len(self.multiframe.microscopes) != len(self.data):
            print("Warning! Data size does not match multiframe model.")
        
        self.recon_type = recon_type
        self.recon_dx = recon_dx
        self.recon_eps = recon_eps
        self.recon_max_iter = recon_max_iter
        self.recon_start = recon_start

        self.recon_est_history = []
        self.recon_norm_history = []
        self.estimated_fluorophore = None
        
    def evaluate(self):
        # Initial estimate
        estimate = self.recon_start
        self.recon_est_history.append(estimate)

        # Iterative procedures
        if self.recon_type == 'Fisher':
            # Fisher scoring algorithm
            iteration = 0
            while iteration < self.recon_max_iter:
                inv_fi = self.noise_model.calc_inv_fi(estimate, 2*[self.recon_dx])
                score = self.noise_model.score(estimate, self.data, 2*[self.recon_dx])
                estimate = estimate + inv_fi.dot(score)
                score_norm = np.linalg.norm(score)
                self.recon_est_history.append(estimate)
                self.recon_norm_history.append(score_norm)
                iteration += 1
                if score_norm < self.recon_eps:
                    break
            self.estimated_fluorophore = fluorophore.Fluorophore(*estimate)
        else:
            # Default reconstruction
            import autograd.numpy as np
            from pymanopt.manifolds import Sphere
            from pymanopt import Problem
            from pymanopt.solvers import ParticleSwarm, NelderMead
            manifold = Sphere(3)
            def cost(X, data=self.data):
                ll = self.multiframe.noise_model.loglikelihood(util.xyz2tp(*X), data)
                return -ll
            problem = Problem(manifold=manifold, cost=cost, verbosity=0)
            start_pts = [np.array(util.tp2xyz(*x)) for x in util.sphere_profile(20)]
            solver = ParticleSwarm(maxcostevals=200)
            Xopt = solver.solve(problem, x=start_pts)
            self.estimated_fluorophore = fluorophore.Fluorophore(*util.xyz2tp(*Xopt))
    
    def plot(self, filename, truth=None, n_pts=1e3):
        # Make axes
        n_cols = 2
        fig = plt.figure()
        fig, axs = plt.subplots(1, n_cols, figsize=(5*n_cols, 5))
        ax_dir = axs[0]
        axs_param = axs[2:-1]
        ax_score = axs[-1]

        estimates = np.array(self.estimates)

        sphere_string = """
        // Sphere
        draw(unitsphere, surfacepen=material(diffusepen=white+opacity(0.1), emissivepen=grey, specularpen=white));
        dotfactor = 7;
        dot(X); 
        dot(Y); 
        circle(0, pi/2, false, (0, 0, 0));
        """
        path_start = 'draw('
        path_middle = ''
        for estimate in estimates[:-1]:
            path_middle += 'expi(theta, phi)--'.replace('theta', str(estimate[0])).replace('phi', str(estimate[1]))
        path_middle += 'expi(theta, phi)'.replace('theta', str(estimates[-1][0])).replace('phi', str(estimates[-1][1]))            
        path_end = ');'
        path_string = path_start + path_middle + path_end
        dot_string = ''
        dot_string += 'dot(expi(theta, phi), black);'.replace('theta', str(estimates[0][0])).replace('phi', str(estimates[0][1]))
        dot_string += 'dot(expi(theta, phi), red);'.replace('theta', str(estimates[-1][0])).replace('phi', str(estimates[-1][1]))
        if truth is not None:
            dot_string += 'dot(expi(theta, phi), green);'.replace('theta', str(truth[0])).replace('phi', str(truth[1]))
        util.draw_scene(dot_string + sphere_string + path_string, my_ax=ax_dir, dpi=300)

        # Plot score_norm
        ax_score.set_yscale('log')
        ax_score.set_xlabel('Iteration Number')
        ax_score.set_ylabel('$||V(\\vec{\\theta})||_2$')
        ax_score.plot(self.score_norm, '-k')

        fig.savefig(filename, dpi=300)
