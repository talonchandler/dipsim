from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

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
                                      n_pts=n_pts_sphere, max_photons=[dose1, dose2], n_samp=1.33)

truth_dirs = util.fibonacci_sphere(1000)
angle_err = []
for truth in truth_dirs:
    data = exp.calc_total_intensities(truth)
    estimate, rh = exp.reconstruct(data, start=np.array([0, 0]), recon_type='Manifold')
    error = np.rad2deg(util.axis_angle(truth, estimate))
    angle_err.append(error)
    print(error, truth, estimate)
    #rh.plot('recon-history3.pdf', truth=truth)

util.plot_sphere('angle-error.pdf', directions=truth_dirs, data=angle_err)
# Plot likelihoood
data = exp.calc_total_intensities((0.01, 0.01))
directions = util.fibonacci_sphere(20000)
I = np.apply_along_axis(exp.noise_model.loglikelihood, 1, directions, data=data)
util.plot_sphere('likelihood.pdf', directions=directions, data=I)

# Draw scene to test
# print('here')
# exp.calc_estimation_stats()

# # Make scene string
# print('here')
# f, ax0 = plt.subplots(1, 1)
# util.draw_scene(exp.scene_string(), my_ax=ax0, dpi=200)
# f.savefig('scene.pdf', dpi=200)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
