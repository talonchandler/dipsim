from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec

# Main input parameters
col_labels = ['Geometry\n (NA = 0.6, $\\beta=80{}^{\circ}$)', 'Uncertainty Ellipses', r'$\sigma_{\Omega}$ [sr]', 'Median$\{\sigma_{\Omega}\}$ [sr]', 'MAD$\{\sigma_{\Omega}\}$ [sr]', '', '']
fig_labels = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)']
n_pts = 5000 #Points on sphere
n_pts_sphere = 50000 # Points on sphere
n_grid_pts = 21
n_line_pts = 50
n_rows, n_cols = 1, len(col_labels)
inch_fig = 5
dpi = 300

# Setup figure and axes
fig = plt.figure(figsize=(2.2*inch_fig, 3*inch_fig))
gs0 = gridspec.GridSpec(3, 1, wspace=0, hspace=0.2, height_ratios=[0.9,1,1])
gs_up = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs0[0], width_ratios=[1, 1, 1, 0.06], wspace=0.1)
gs_middle = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[1], width_ratios=[1, 1], wspace=0.4)
gs_down = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs0[2], width_ratios=[1, 1], wspace=0.4)
gs_middle_left = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_middle[0], width_ratios=[1, 0.05], wspace=0.1)
gs_middle_right = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_middle[1], width_ratios=[1, 0.05], wspace=0.1)
gs_down_left = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_down[0], width_ratios=[1, 0.05], wspace=0.1)
gs_down_right = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_down[1], width_ratios=[1, 0.05], wspace=0.1)

ax0 = plt.subplot(gs_up[0])
ax1 = plt.subplot(gs_up[1])
ax2 = plt.subplot(gs_up[2])
cax2 = plt.subplot(gs_up[3])
ax3 = plt.subplot(gs_middle_left[0])
cax3 = plt.subplot(gs_middle_left[1])
ax4 = plt.subplot(gs_middle_right[0])
cax4 = plt.subplot(gs_middle_right[1])
ax5 = plt.subplot(gs_down_left[0])
cax5 = plt.subplot(gs_down_left[1]); cax5.axis('off');
ax6 = plt.subplot(gs_down_right[0])
cax6 = plt.subplot(gs_down_right[1]); cax6.axis('off');

for ax, col_label, fig_label  in zip([ax0, ax1, ax2, ax3, ax4, ax5, ax6], col_labels, fig_labels):
    ax.annotate(col_label, xy=(0,0), xytext=(0.5, 1.06), textcoords='axes fraction',
                va='center', ha='center', fontsize=14, annotation_clip=False)
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 1.06), textcoords='axes fraction',
                va='center', ha='center', fontsize=14, annotation_clip=False)
    
for ax in [ax0, ax1, ax2, ax3, ax4, ax5, ax6]:
    ax.tick_params(axis='both', labelsize=14)
for cax in [cax2, cax3, cax4]:
    cax.tick_params(axis='both', labelsize=14)

# Calculate a list of points to sample in region
n = 1.33
NA_max = n*np.sin(np.pi/4)
NA = np.linspace(0, n, 1000)
lens_bound = np.rad2deg(2*np.arcsin(NA/n))
cover_bound = np.rad2deg(np.pi - 2*np.arcsin(NA/n))

pts = np.mgrid[n/n_grid_pts/2:n:n_grid_pts*1j,0:180:n_grid_pts*1j].reshape((2, n_grid_pts**2)).T.tolist()

def is_feasible(pt):
    if pt[1] < np.rad2deg(np.pi - 2*np.arcsin(pt[0]/n)) + 20 and pt[1] > np.rad2deg(2*np.arcsin(pt[0]/n)) - 20 and pt[0] < NA_max + 0.1:
        return True
    else:
        return False

pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate med and mad for each point
def calc_stats(param):
    na = param[0]
    angle = param[1]
    exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.deg2rad(angle/2), -np.deg2rad(angle/2)], det_thetas=[-np.deg2rad(angle/2), np.deg2rad(angle/2)],
                                      ill_nas=2*[na], det_nas=2*[na],
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
def plot_2d_regions(ax, cax, pts, data, special_pt=(-1,-1),
                    line_pt0=None, line_pt1=None):
    ax.plot(NA, lens_bound, 'k-', zorder=11)
    ax.plot(NA, cover_bound, 'k-', zorder=11)

    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator
    def degrees(x, pos):
        return str(int(x)) + '${}^{\circ}$'
    ax.yaxis.set_major_locator(FixedLocator([0, 45, 90, 135, 180]))
    ax.yaxis.set_major_formatter(FuncFormatter(degrees))

    from matplotlib.ticker import FuncFormatter, FixedLocator
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0, 1.33])
    ax.set_xticklabels(['0', '0.25', '0.5', '0.75', '1.0', '1.33'])

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, 'NA', (0.5, -0.12), fontsize=14)
    my_annotate(ax, '$\\beta$, Angle Between Objectives', (-0.25, 0.5), fontsize=14, rotation=90)
    my_annotate(ax, 'Objectives collide\nwith cover slip', (0.65, 0.85), fontsize=14)
    my_annotate(ax, 'Objectives collide\nwith each other', (0.65, 0.15), fontsize=14)
    my_annotate(ax, 'Feasible', (0.3, 0.5), fontsize=14)

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
    
    ax.plot([line_pt0[0], line_pt1[0]], [line_pt0[1], line_pt1[1]], '-', color='darkmagenta', lw=3, zorder=1)
    ax.plot(special_pt[0], special_pt[1], 'kx', markersize=5)


    # Plot patches
    width = n/(n_grid_pts)
    for i, (pt, c) in enumerate(zip(pts_list, colors)):
        if pt[1] == 0:
            height = 180/14.5
        if pt[1] == 0:
            height = 180/(n_grid_pts-0.5)
        ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
    fig.colorbar(sc, cax=cax, orientation='vertical')

    # Mask around lines
    ax.fill_between(NA, lens_bound, 0, color='white', zorder=2)
    ax.fill_between(NA, cover_bound, 180, color='white', zorder=2)

    ax.set(xlim=[0, 1.33], ylim=[0, 180])

# Plot 1D region
def plot_1d_regions(ax, pts, data, special_pt=(-1,-1), y_pos=None, y_lim=None, xtitle=None):
    # Set y ticks
    from matplotlib.ticker import FuncFormatter, FixedLocator
    def degrees(x, pos):
        return str(int(x)) + '${}^{\circ}$'
    ax.xaxis.set_major_locator(FixedLocator([53, 90, 135, 127]))
    ax.xaxis.set_major_formatter(FuncFormatter(degrees))

    from matplotlib.ticker import FuncFormatter, FixedLocator
    ax.set_yticks(y_pos)
    ax.set_yticklabels(["{:.1e}".format(x).replace('e-0', 'e-') for x in y_pos])

    # Annotation
    def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
        ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                    va='center', ha='center', fontsize=fontsize,
                    annotation_clip=False, rotation=rotation, zorder=13)

    my_annotate(ax, '$\\beta$, Angle Between Objectives', (0.5, -0.12), fontsize=14)
    my_annotate(ax, xtitle, (-0.25, 0.5), fontsize=14, rotation=90)

    ax.set(xlim=[53, 127], ylim=y_lim)
    ax.plot(pts, data, '-', color='darkmagenta', lw=3, zorder=1)
    
# Plot first two columns
angle = 80
na = 0.6
exp = multiframe.MultiFrameMicroscope(ill_thetas=[np.deg2rad(angle/2), -np.deg2rad(angle/2)], det_thetas=[-np.deg2rad(angle/2), np.deg2rad(angle/2)],
                                  ill_nas=2*[na], det_nas=2*[na],
                                  ill_types=2*['wide'], det_types=2*['lens'],
                                  colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                  n_pts=n_pts_sphere, max_photons=500, n_samp=1.33)
exp.calc_estimation_stats()

# Make scene string
scene_string = exp.scene_string()
line_string = "draw(O--expi(theta, 0));\n"
line_string = line_string.replace('theta', str(np.deg2rad(angle/2)))
scene_string += line_string
line_string = "draw(O--expi(theta, 0));\n"
line_string = line_string.replace('theta', str(np.deg2rad(-angle/2)))
scene_string += line_string
arc_string = 'draw(L=Label("$\\beta$", align=N), arc(O, 0.1*expi(-theta, 0), 0.1*expi(theta, 0), normal=Y));\n'
arc_string = arc_string.replace('theta', str(np.deg2rad(angle/2)))
scene_string += arc_string

util.draw_scene(scene_string, my_ax=ax0, dpi=dpi)
util.draw_scene(exp.ellipse_string(n_pts=250), my_ax=ax1, dpi=dpi)
util.plot_sphere(directions=exp.directions, data=exp.sa_uncert,
                 color_norm='log', linthresh=1e-4,
                 color_min=1e-4, color_max=1e1,
                 my_ax=ax2, my_cax=cax2)

# Find profile points
line_na = 0.6
min_beta = np.rad2deg(2*np.arcsin(line_na/n))
max_beta = 180 - np.rad2deg(2*np.arcsin(line_na/n))

# Plots last two columns
plot_2d_regions(ax3, cax3, pts, med, special_pt=(na, angle), line_pt0=(line_na, min_beta), line_pt1=(line_na, max_beta))
plot_2d_regions(ax4, cax4, pts, mad, special_pt=(na, angle), line_pt0=(line_na, min_beta), line_pt1=(line_na, max_beta))

# Calculate and plot profile
line_beta = np.linspace(min_beta, max_beta, n_line_pts)
line_na = 0.6*np.ones(line_beta.shape)
line_pts = np.vstack([line_na, line_beta])
line_med = []
line_mad = []
for i, pt in enumerate(line_pts.T):
    print('Calculating microscope '+str(i+1)+'/'+str(line_pts.shape[1]))
    x = calc_stats(pt)
    line_med.append(x[0])
    line_mad.append(x[1])

plot_1d_regions(ax5, line_beta, line_med, special_pt=angle, y_pos=[4.5e-3, 5e-3, 5.5e-3], y_lim=[4.4e-3, 5.6e-3], xtitle='Median$\{\sigma_{\Omega}\}$ [sr]')
plot_1d_regions(ax6, line_beta, line_mad, special_pt=angle, y_pos=[1e-3, 1.5e-3, 2e-3], y_lim=[8e-4, 2e-3], xtitle='MAD$\{\sigma_{\Omega}\}$ [sr]')

# Label axes and save
print('Saving final figure.')    
fig.savefig('../paper/symmetric-widefield.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
