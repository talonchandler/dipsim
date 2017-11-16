import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from dipsim import util, multiframe, fluorophore
from pymanopt.manifolds import Product, Sphere, Euclidean
from pymanopt import Problem
from pymanopt.solvers import ParticleSwarm

class Reconstruction():
    """A reconstruction is specified by a MuliFrameMicroscope, data, and a 
    reconstruction options. 
    """
    def __init__(self, multiframe=multiframe.MultiFrameMicroscope(),
                 data=None, recon_type='tp', recon_dx=1e-6, recon_eps=1e-3,
                 recon_max_iter=1e2):
        
        self.multiframe = multiframe
        self.data = data
        if len(self.multiframe.microscopes) != len(self.data):
            print("Warning! Data size does not match multiframe model.")
        
        self.recon_type = recon_type
        self.recon_dx = recon_dx
        self.recon_eps = recon_eps
        self.recon_max_iter = recon_max_iter

        self.recon_est_history = []
        self.recon_norm_history = []
        self.estimated_fluorophore = None
        
    def evaluate(self):
        # Perform reconstruction
        if self.recon_type == 'tp': # Estimate theta and phi
            manifold = Sphere(3)
            def cost(X, data=self.data):
                ll = self.multiframe.noise_model.loglikelihood(util.xyz2tp(*X), data)
                return -ll
            problem = Problem(manifold=manifold, cost=cost, verbosity=0)
            start_pts = [np.array(util.tp2xyz(*x)) for x in util.sphere_profile(20)]
            solver = ParticleSwarm(maxcostevals=200)
            Xopt = solver.solve(problem, x=start_pts)
            self.estimated_fluorophore = fluorophore.Fluorophore(*util.xyz2tp(*Xopt))
        elif self.recon_type == 'tpc': # Estimate theta, phi, constant
            # Create manifold and cost function
            manifold = Product((Sphere(3), Euclidean(1)))
            def cost(X, data=self.data):
                estimate = np.hstack([util.xyz2tp(*X[0]), X[1]])
                ll = self.multiframe.noise_model.loglikelihood(estimate, data)
                return -ll
            
            problem = Problem(manifold=manifold, cost=cost, verbosity=0)

            # Generate start_pts and format            
            xyz_start_pts = 3*[np.array(util.tp2xyz(*x)) for x in util.sphere_profile(10)]
            c_start_pts = np.expand_dims(np.hstack((10*[0.1], 10*[2], 10*[10])), axis=1)
            start_pts = np.hstack((xyz_start_pts, c_start_pts))
            pts = []
            for start_pt in start_pts:
                pts.append([np.array(start_pt[0:3]), np.array(start_pt[3:5])])

            # Solve
            solver = ParticleSwarm(maxcostevals=500)
            Xopt = solver.solve(problem, x=pts)

            self.estimated_fluorophore = fluorophore.Fluorophore(*np.hstack([util.xyz2tp(*Xopt[0]), Xopt[1]]).flatten())

        elif self.recon_type == 'tpck': # Estimate theta, phi, constant, kappa
            # Create manifold and cost function
            manifold = Product((Sphere(3), Euclidean(2)))
            def cost(X, data=self.data):
                estimate = np.array([util.xyz2tp(*X[0]), X[1]]).flatten() # Reshape data for loglikelihood function
                ll = self.multiframe.noise_model.loglikelihood(estimate, data)
                print(estimate, ll)
                return -ll
            
            problem = Problem(manifold=manifold, cost=cost, verbosity=0)

            # Generate start_pts and format            
            xyz_start_pts = 3*[np.array(util.tp2xyz(*x)) for x in util.sphere_profile(10)]
            k_start_pts = np.expand_dims(np.hstack((10*[-100], 10*[0], 10*[100])), axis=1)
            c_start_pts = np.expand_dims(np.hstack((10*[0.1], 10*[1], 10*[10])), axis=1)
            start_pts = np.hstack((xyz_start_pts, c_start_pts, k_start_pts))
            pts = []
            for start_pt in start_pts:
                pts.append([np.array(start_pt[0:3]), np.array(start_pt[3:5])])

            # Solve
            solver = ParticleSwarm(maxcostevals=200)
            Xopt = solver.solve(problem, x=pts)

            self.estimated_fluorophore = fluorophore.Fluorophore(*np.array([util.xyz2tp(*Xopt[0]), Xopt[1]]).flatten())
    
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
