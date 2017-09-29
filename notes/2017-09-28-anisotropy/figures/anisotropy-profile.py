from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Main input parameters
n_pts = 200
kappas = [-3, 0, 3, None]
n_cols = 3
n_rows = 1
inch_fig = 5
dpi = 300
col_labels = ['$I_\perp$', '$I_\parallel$', '$r = \\frac{I_\parallel - I_\perp}{I_\parallel + 2I_\perp}$']

# Generate axes
size = (inch_fig*n_cols, 0.7*inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.3, hspace=-0.1)
#caxs = util.generate_caxs(axs)

# Compute and plot on axes
for i, kappa in enumerate(kappas):
    print('Computing microscope: ' + str(i))

    # y-polarized input from x optical axis
    ill = illuminator.Illuminator(illum_type='sheet', theta_optical_axis=np.pi/2, phi_pol=np.pi/2)
    # perp: x-polarized detection
    # par: y-polarized detection
    det_perp = detector.Detector(det_type='polarized', na=0.006, phi_pol=0)
    det_par = detector.Detector(det_type='polarized', na=0.006, phi_pol=np.pi/2)
    m_perp = microscope.Microscope(illuminator=ill, detector=det_perp, max_photons=1e6)
    m_par = microscope.Microscope(illuminator=ill, detector=det_par, max_photons=1e6)
    
    directions = util.sphere_profile(n_pts)
    I_perp = np.apply_along_axis(m_perp.calc_intensity, 1, directions, kappa=kappa)
    I_par = np.apply_along_axis(m_par.calc_intensity, 1, directions, kappa=kappa)
    r = (I_par - I_perp)/(I_par + 2*I_perp)
    
    phi = directions[:,1]
    if kappa == None:
        k_string = '$\infty$'
    else:
        k_string = str(kappa)
    axs[0].plot(np.rad2deg(phi), I_perp, '-')
    axs[1].plot(np.rad2deg(phi), I_par, '-', label='$\kappa = \ $'+k_string)
    axs[2].plot(np.rad2deg(phi), r, '-')

axs[0].set_ylim([0, 8])    
axs[1].set_ylim([0, 8])
axs[2].set_ylim([-0.5, 1.0])
axs[0].set_ylabel('$I_\perp$')
axs[1].set_ylabel('$I_\parallel$')
axs[2].set_ylabel('$r = \\frac{I_\parallel - I_\perp}{I_\parallel + 2I_\perp}$')
axs[1].legend(loc='upper center')



for ax in axs:
    ax.set_xlim([-180, 180])
    ax.set_xlabel('$\phi$')
    from matplotlib.ticker import FuncFormatter, FixedLocator
    def degrees(x, pos):
        return str(int(x)) + '${}^{\circ}$'
    ax.xaxis.set_major_locator(FixedLocator([-180, -90, 0, 90, 180]))
    ax.xaxis.set_major_formatter(FuncFormatter(degrees))
    
    
# Label axes and save
#util.label_rows_and_cols(axs, row_labels, col_labels)
print('Saving final figure.')    
fig.savefig('anisotropy-profile.pdf', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
