from dipsim import multiframe, util, detector, illuminator, microscope, util
import dipsim.fluorophore as flu
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 1000
kappas = [-3, 0, 3, None]
n_cols = len(kappas)
n_rows = 3
inch_fig = 5
dpi = 300

col_labels = ['$\kappa = -3$', '$\kappa = 0$', '$\kappa = 3$', '$\kappa = \infty$']
row_labels = ['$I_\perp$', '$I_\parallel$', '$r = \\frac{I_\parallel - I_\perp}{I_\parallel + 2I_\perp}$']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, kappa in enumerate(kappas):
    print('Computing microscope: ' + str(i))

    # y-polarized input from x optical axis
    ill = illuminator.Illuminator(illum_type='sheet', theta_optical_axis=np.pi/2, phi_pol=np.pi/2)
    # perp: x-polarized detection
    # par: y-polarized detection
    det_perp = detector.Detector(det_type='polarized', na=0.01, phi_pol=0)
    det_par = detector.Detector(det_type='polarized', na=0.01, phi_pol=np.pi/2)
    m_perp = microscope.Microscope(illuminator=ill, detector=det_perp, max_photons=1e6)
    m_par = microscope.Microscope(illuminator=ill, detector=det_par, max_photons=1e6)
    test_perp = m_perp.calc_intensity((0,0), kappa=0)
    test_par = m_par.calc_intensity((0,0), kappa=0)    
    directions = util.fibonacci_sphere(n_pts)
    I_perp = np.apply_along_axis(m_perp.calc_intensity, 1, directions, kappa=kappa)
    I_par = np.apply_along_axis(m_par.calc_intensity, 1, directions, kappa=kappa)
    
    r = (I_par - I_perp)/(I_par + 2*I_perp)
        
    util.plot_sphere('', directions=directions, data=I_perp, my_ax=axs[0,i], my_cax=caxs[0,i],
                     color_norm='log', color_min=1e-1, color_max=3e1)
    util.plot_sphere('', directions=directions, data=I_par, my_ax=axs[1,i], my_cax=caxs[1,i],
                     color_norm='log', color_min=1e-1, color_max=3e1)                     
    util.plot_sphere('', directions=directions, data=r, my_ax=axs[2,i], my_cax=caxs[2,i],
                     color_min=-0.5, color_max=1.0)

# Label axes and save
util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('anisotropy-kappa.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
