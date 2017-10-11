import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from dipsim import util

class ReconHistory:
    def __init__(self):
        self.truth = None
        self.estimates = []        
        self.score_norm = []
        self.eps = 0

    def plot(self, filename, iter_plot_n=3, truth=None, n_pts=1e3):
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
