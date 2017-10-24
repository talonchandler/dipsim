from dipsim import multiframe, util, fluorophore, reconstruction
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Setup k and c sweep
kappas = [10 , 100, np.inf]
cs = [0.1, 1, 2]
col_labels = ['$\kappa$ = ' + str(x) for x in kappas]
row_labels = ['c = ' + str(x) for x in cs]

# Generate axes
inch_fig = 5
n_rows = len(cs)
n_cols = len(kappas)
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0)
if len(col_labels) == 1:
    axs = np.expand_dims(axs, 1)
caxs = util.generate_caxs(axs)

# Specify microscope geometry
n = 1.33
alpha1 = 60
na1 = 1.1
na2 = 0.71
theta1 = np.deg2rad(45-12)
theta2 = -np.deg2rad(45+12)
dose_ratio = 1
total_photons = 10000
dose2 = total_photons/(1 + dose_ratio)
dose1 = total_photons - dose2

exp = multiframe.MultiFrameMicroscope(ill_thetas=[theta1, theta2], det_thetas=[theta2, theta1],
                                      ill_nas=[na1, na2], det_nas=[na2, na1],
                                      ill_types=2*['sheet'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      max_photons=[dose1, dose2], n_samp=1.33)

#exp.precompute_efficiencies(n_pts=1e2)

# Plot likelihoood
for i, kappa in enumerate(kappas):
    for j, c in enumerate(cs):
        print(i, j)
        fluo = fluorophore.Fluorophore(theta=0.01, phi=0.01)
        data = exp.calc_intensities(fluo)
        n_pts = 50
        directions = util.fibonacci_sphere(n_pts)
        input_points = np.hstack([directions, kappa*np.ones([n_pts, 1]), c*np.ones([n_pts, 1])])
        likelihoods = [exp.noise_model.loglikelihood(x, data=data) for x in input_points]
        util.plot_sphere(directions=directions, data=likelihoods,
                         my_ax=axs[j,i], my_cax=caxs[j,i],
                         color_min=1e1, color_max=1e3)

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels, row_pos=(-0.2, 0.5))
print('Saving final figure.')    
fig.savefig('likelihood.pdf', dpi=300)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
