from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry (NA${}_{\\textrm{ill}}$ = 0.2, NA${}_{\\textrm{det}}$ = 0.8\n 3 illumination polarizations)', 'Uncertainty Ellipse', r'$\sigma_{\Omega}$ [sr]', 'Median$\{\sigma_{\Omega}\}$ [sr]', 'MAD$\{\sigma_{\Omega}\}$ [sr]']
fig_labels = ['a)', 'b)', 'c)', 'd)', 'e)']
n_pts = 1000 # Points on sphere
n_pts_sphere = 20000 # Points on sphere
n_grid_pts = 10
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(2.2*inch_fig, 2*inch_fig))
gs0 = gridspec.GridSpec(2, 1, wspace=0, hspace=0.2, height_ratios=[0.9,1])
gs_up = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[0], width_ratios=[1, 1, 1, 0.06], wspace=0.1)
gs_middle = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], width_ratios=[1, 1], wspace=0.4)
gs_middle_left = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_middle[0], width_ratios=[1, 0.05], wspace=0.1)
gs_middle_right = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_middle[1], width_ratios=[1, 0.05], wspace=0.1)

ax0 = plt.subplot(gs_up[0])
ax1 = plt.subplot(gs_up[1])
ax2 = plt.subplot(gs_up[2])
cax2 = plt.subplot(gs_up[3])
ax3 = plt.subplot(gs_middle_left[0])
cax3 = plt.subplot(gs_middle_left[1])
ax4 = plt.subplot(gs_middle_right[0])
cax4 = plt.subplot(gs_middle_right[1])

for ax, col_label, fig_label  in zip([ax0, ax1, ax2, ax3, ax4], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
    
for ax in [ax0, ax1, ax2, ax3, ax4]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax2, cax3, cax4]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 1.33
NA_max = n
x = np.arange(2, 11, 1)
y, step = np.linspace(0, NA_max, num=n_grid_pts, retstep=True, endpoint=False)
y += step/2
pts = np.array(np.meshgrid(x, y)).reshape(2, len(x)*len(y)).T

def is_feasible(pt):
    return True

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate med and mad for each point
def calc_stats(param):
    n_pol = int(param[0])
    ill_na = param[1]
    det_na = 0.8

    exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                          ill_nas=2*[ill_na], det_nas=2*[det_na],
                                          ill_types=2*['wide'], det_types=2*['lens'],
                                          colors=['(1,0,0)', '(0,0,1)'], n_frames=n_pol,
                                          n_pts=n_pts, max_photons=2000/n_pol, n_samp=1.33)
    
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
    
    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, '\# of illumination polarizations', (0.5, -0.13), fontsize=14)
    my_annotate(ax, 'NA${}_{\\textrm{ill}}$', (-0.13, 0.5), fontsize=14, rotation=90)

    # Calculate colors
    color_map='coolwarm'
    color_norm='log'
    color_min=1e-4
    color_max=1e1
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
    width = 1
    height = NA_max/len(y)
    for i, (pt, c) in enumerate(zip(pts_list, colors)):
        ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
    fig.colorbar(sc, cax=cax, orientation='vertical')

    ax.set(xlim=[1.5, 10], ylim=[0, 1.33])

# Plot first two columns
n_pol = 3
ill_na = 0.2
det_na = 0.8
exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.pi/4, -np.pi/4], det_thetas=[-np.pi/4, np.pi/4],
                                      ill_nas=2*[ill_na], det_nas=2*[det_na],
                                      ill_types=2*['wide'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=n_pol,
                                      n_pts=n_pts, max_photons=2000/n_pol, n_samp=1.33)

exp.calc_estimation_stats()

# Make scene string

util.draw_scene(exp.scene_string(), my_ax=ax0, dpi=dpi)
util.draw_scene(exp.ellipse_string(n_pts=250), my_ax=ax1, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_max=1e1, color_min=1e-4,
                 my_ax=ax2, my_cax=cax2)
    
# Plots last two columns
plot_2d_regions(ax3, cax3, pts, med, special_pt=(n_pol, ill_na))
plot_2d_regions(ax4, cax4, pts, mad, special_pt=(n_pol, ill_na))
    
# Label axes and save
print('Saving final figure.')    
fig.savefig('ill-polarization-ortho.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
