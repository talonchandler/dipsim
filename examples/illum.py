from dipsim import fluorophore, illuminator, detector, microscope, stats, util
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

from mpl_toolkits.axes_grid1 import make_axes_locatable

n = 200#2000
dpi = 250
vis_px = 2000
bfp_rads = [1, 2, 3, 4]
n_cols = len(bfp_rads)
n_rows = 5
inch_fig = 5

fig, axs = plt.subplots(n_rows, n_cols,
                        figsize=(inch_fig*n_cols, inch_fig*n_rows))

caxs = []
for ax in axs.flatten():
    divider = make_axes_locatable(ax)
    caxs.append(divider.append_axes("right", size="5%", pad=0.15))
caxs = np.array(caxs).reshape(axs.shape)
    
if len(bfp_rads) == 1:
    axs = np.expand_dims(axs, 1)

plt.subplots_adjust(wspace=0.2, hspace=0)

for i, bfp_rad in enumerate(bfp_rads):
    f = 10
    bfp_n = 8

    illx = illuminator.Illuminator(illum_type='kohler',
                                   optical_axis=np.array([0., 0., 1.]),
                                   f=f,
                                   bfp_rad=bfp_rad,
                                   bfp_pol=np.array([1., 0., 0.]),
                                   bfp_n=bfp_n)

    det = detector.Detector(optical_axis=np.array([0, 0, 1]),
                            na=1.5,
                            n=1.5)

    m = microscope.Microscope(illuminator=illx, detector=det)

    m.plot_intensities_from_single_fluorophore('intensity'+str(bfp_rad)+'.png',
                                                n=n, interact=False, color_norm='linear', my_ax=axs[1][i], my_cax=caxs[1][i], dpi=dpi, vis_px=vis_px)
    
    m.draw_scene('scene'+str(bfp_rad)+'.png', interact=False, my_ax=axs[0][i], dpi=dpi, vis_px=vis_px)
    caxs[0][i].axis('off')

row_labels = ['Scene', 'Emitted Power', r'$\sigma_{\phi}$', r'$\sigma_{\theta}$', r'$\sigma_{\Omega} = \sigma_{\phi}\sigma_{\theta}\sin\theta$']
for i, label in enumerate(row_labels):
    axs[i][0].annotate(label, xy=(0,0), xytext=(-0.1, 0.5), textcoords='axes fraction',
                       va='center', ha='center', rotation=90, fontsize=18)

col_labels = [r'$N=1$', r'$N=2$', r'$N=3$', r'$N=4$']
for i, label in enumerate(col_labels):
    axs[0][i].annotate(label, xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction',
                       va='center', ha='center', fontsize=18)

fig.savefig('compare_N_pol.png', dpi=dpi)
print('Total time: '+str(np.round(time.time() - start, 2)))
