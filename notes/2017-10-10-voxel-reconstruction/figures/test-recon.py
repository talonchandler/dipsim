from dipsim import multiframe, util, fluorophore, reconstruction
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Specify microscope geometry
n_pts_sphere = 500 # Points on sphere

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


# Test reconstructions
truth_dirs = util.fibonacci_sphere(250)
angle_err = []
for truth in truth_dirs:
    fluo = fluorophore.Fluorophore(theta=truth[0], phi=truth[1])
    recon = reconstruction.Reconstruction(multiframe=exp, data=exp.calc_intensities(fluo))
    recon.evaluate()
    error = fluo - recon.estimated_fluorophore
    angle_err.append(error['angle_diff'])
    print(error['angle_diff'], truth)

util.plot_sphere('angle-error.pdf', directions=truth_dirs, data=angle_err)

# Plot likelihoood
fluo = fluorophore.Fluorophore(theta=0.01, phi=0.01)
data = exp.calc_intensities(fluo)
directions = util.fibonacci_sphere(20000)
likelihoods = [exp.noise_model.loglikelihood(x, data=data) for x in directions]
util.plot_sphere('likelihood.pdf', directions=directions, data=likelihoods)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
