from dipsim import illuminator, util, microscope, detector
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')

# Main input parameters
thetas = [x*np.pi/4 for x in range(4)]
bfp_rads = np.hstack((np.arange(0.01, 0.1, 0.01), np.arange(0.2, 11, 0.5)))
n_cols = len(thetas)
n_rows = 3
inch_fig = 5
vis_px = 2000
dpi = 300
row_labels = ['Scene']

# Generate axes
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.2, hspace=0.4)
if len(thetas) == 1:
    axs = np.expand_dims(axs, 1)

# Compute and plot on axes
for i, theta in enumerate(thetas):
    print('Computing microscope: ' + str(theta))
    ill_basis = []
    for bfp_rad in bfp_rads:
        bfp_pol = np.array([np.cos(theta), np.sin(theta), 0])
        f = 10
        ill = illuminator.Illuminator(illum_type='kohler',
                                      optical_axis=np.array([0., 0., 1.]),
                                      f=10, bfp_rad=bfp_rad,
                                      bfp_pol=bfp_pol,
                                      bfp_n=256) #>128 for good
        ill_basis.append(ill.illum_basis)

    ill_basis = np.array(ill_basis)

    for ax_i in [1, 2]:
        for j in range(ill_basis.shape[1]):
            axs[ax_i,i].plot(bfp_rads/f, ill_basis[:,j], '-', lw=1.0+1.0*((j+1)%2), alpha=0.8)
        labels = [r'$|A_x|^2$', r'$|A_y|^2$', r'$|A_z|^2$', r'2Re\{$A_x^*A_y$\}', r'2Re\{$A_x^*A_z$\}', r'2Re\{$A_y^*A_z$\}']
        axs[ax_i,i].legend(labels, bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=3, mode="expand", borderaxespad=0.)
        axs[ax_i,i].set_xlabel(r'\rho/f')
        axs[ax_i,i].set_xlim([0,1.0])
        if ax_i == 1:
            axs[ax_i,i].set_yscale('linear')
            axs[ax_i,i].set_ylim([-1.0,1.0])
        else:
            axs[ax_i,i].set_yscale('log')
            axs[ax_i,i].set_ylim([1e-5,1.0])
    
    det = detector.Detector(optical_axis=np.array([0, 0, 1]), na=1.5, n=1.5)
    ill.bfp_rad = 3
    m = microscope.Microscope(illuminator=ill, detector=det, max_photons=1)
    m.draw_scene(my_ax=axs[0,i], dpi=dpi, vis_px=vis_px)

# Label axes and save
util.label_rows_and_cols(axs, row_labels)
print('Saving final figure.')    
fig.savefig('illum_basis.png', dpi=dpi)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
