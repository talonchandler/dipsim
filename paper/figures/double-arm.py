from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry\n(NA${}_{\\textrm{ill}}$ = 0.6, NA${}_{\\textrm{det}}$ = 0.8)', r'$\sigma_{\Omega}$ [sr]', 'Median$\{\sigma_{\Omega}\}$ [sr]', 'MAD$\{\sigma_{\Omega}\}$ [sr]']
fig_labels = ['a)', 'b)', 'c)', 'd)']
n_pts = 500 # Points on sphere
n_pts_sphere = 50000 # Points on sphere
n_grid_pts = 25
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(2.2*inch_fig, 2*inch_fig))
gs0 = gridspec.GridSpec(2, 2, wspace=0.4, hspace=0.1)
gs00 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,0], width_ratios=[1, 0.05], wspace=0.1)
gs10 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,0], width_ratios=[1, 0.05], wspace=0.1)
gs01 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[0,1], width_ratios=[1, 0.05], wspace=0.1)
gs11 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1,1], width_ratios=[1, 0.05], wspace=0.1)
ax0 = plt.subplot(gs00[0])
ax1 = plt.subplot(gs01[0])
ax2 = plt.subplot(gs10[0])
ax3 = plt.subplot(gs11[0])
cax0 = plt.subplot(gs00[1]); cax0.axis('off');
cax1 = plt.subplot(gs01[1])
cax2 = plt.subplot(gs10[1])
cax3 = plt.subplot(gs11[1])

for ax, col_label, fig_label  in zip([ax0, ax1, ax2, ax3], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    
for ax in [ax0, ax1, ax2, ax3]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax1, cax2, cax3]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 0.94
NA_ill = np.linspace(0, 0.94, num=1000)
NA_det = NA_ill
pts = np.mgrid[n/n_grid_pts/2:n-n/n_grid_pts/2:n_grid_pts*1j,n/n_grid_pts/2:n-n/n_grid_pts/2:n_grid_pts*1j].reshape((2, n_grid_pts**2)).T.tolist()

def is_feasible(pt):
    if pt[0] < pt[1] + 0.2:
        return True
    else:
        return False

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate med and mad for each point
def calc_stats(param):
    na_ill = param[0]
    na_det = param[1]
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                      ill_nas=2*[na_ill], det_nas=2*[na_det],
                                      ill_types=2*['wide'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      n_pts=n_pts, max_photons=500, n_samp=1.33)

    exp.calc_estimation_stats()
    data = exp.sa_uncert
    med = np.median(data)
    return med, np.median(np.abs(data - med))

med = []
mad = []
for i, pt in enumerate(pts.T):
    print('Calculating microscope '+str(i+1)+'/'+str(pts.shape[1]))
    x = calc_stats(pt)
    med.append(x[0])
    mad.append(x[1])

# Plot 2D regions
def plot_2d_regions(ax, cax, pts, data, special_pt=(-1,-1)):
    ax.plot(NA_ill, NA_det, 'k-', zorder=11)

    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator
    ax.set_yticks([0, 0.25, 0.5, 0.75, 0.94])
    ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '0.94'])
    ax.set_xticks([0, 0.25, 0.5, 0.75, 0.94])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '0.94'])

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, 'NA${}_{\\textrm{ill}}$', (0.5, -0.15), fontsize=14)
    my_annotate(ax, 'NA${}_{\\textrm{det}}$', (-0.18, 0.5), fontsize=14, rotation=90)

    # Calculate colors
    color_map='coolwarm'
    color_norm='log'
    color_min=1e-4
    color_max=1e-1
    if color_norm == 'linear':
        norm = matplotlib.colors.Normalize(vmin=color_min, vmax=color_max)
    elif color_norm == 'log':
        norm = matplotlib.colors.LogNorm(vmin=color_min, vmax=color_max)
    elif color_norm == 'linlog':
        norm = matplotlib.colors.SymLogNorm(linthresh=linthresh, vmin=-color_max, vmax=color_max)
    elif color_norm == 'power':
        norm = matplotlib.colors.PowerNorm(gamma=gamma, vmin=data.min(), vmax=data.max())

    norm_data = norm(data).data
    norm_data2 = np.expand_dims(norm_data, 1)    
    cmap = matplotlib.cm.get_cmap(color_map)
    colors = np.apply_along_axis(cmap, 1, norm_data2)

    # Plot scatter for colorbar
    sc = ax.scatter(pts[0,:], pts[1,:], c=data, s=0, cmap=cmap, norm=norm,
                    marker='s', lw=0)

    ax.plot(special_pt[0], special_pt[1], 'kx', markersize=5, zorder=15)

    # Plot patches
    width = n/(n_grid_pts)
    height = width
    for i, (pt, c) in enumerate(zip(pts_list, colors)):
        ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
    fig.colorbar(sc, cax=cax, orientation='vertical')

    # Mask around lines
    ax.fill_between(NA_ill, NA_det, 0, color='white', zorder=2)
    ax.set(xlim=[0, 0.94], ylim=[0, 0.94])

# Plot first two columns
na_ill = 0.6
na_det = 0.8
exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                  ill_nas=2*[na_ill], det_nas=2*[na_det],
                                  ill_types=2*['wide'], det_types=2*['lens'],
                                  colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                  n_pts=n_pts_sphere, max_photons=500, n_samp=1.33)

exp.calc_estimation_stats()

# Make scene string
scene_string = exp.scene_string()
util.draw_scene(scene_string, my_ax=ax0, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=None, color_max=None,
                 my_ax=ax1, my_cax=cax1)
    
# Plots last two columns
plot_2d_regions(ax2, cax2, pts, med, special_pt=(na_ill, na_det))
plot_2d_regions(ax3, cax3, pts, mad, special_pt=(na_ill, na_det))
    
# Label axes and save
print('Saving final figure.')    
fig.savefig('../paper/double-arm.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
