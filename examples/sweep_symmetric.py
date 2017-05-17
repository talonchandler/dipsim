from dipsim import multiframe, util
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches as patches

import os; import time; start = time.time(); print('Running...')

# Plotting
n = 1.33
NA_max = n*np.sin(np.pi/4)
NA = np.linspace(0, n, 1000)
lens_bound = np.rad2deg(2*np.arcsin(NA/n))
cover_bound = np.rad2deg(np.pi - 2*np.arcsin(NA/n))
f = plt.figure(figsize=(5,5))
ax = plt.axes([0, 0, 1.0, 1.0]) # x, y, width, height
ax.set_aspect(n/180)
cax = plt.axes([1.05, 0, 0.025, 1.0])
ax.plot(NA, lens_bound, 'k-', zorder=11)
ax.plot(NA, cover_bound, 'k-', zorder=11)
ax.set_xlim([0, 1.33])
ax.set_ylim([0, 180])

# Set y ticks
from matplotlib.ticker import FuncFormatter, FixedLocator
def degrees(x, pos):
    return str(int(x)) + '${}^{\circ}$'
ax.yaxis.set_major_locator(FixedLocator([0, 45, 90, 135, 180]))
ax.yaxis.set_major_formatter(FuncFormatter(degrees))

# Annotation
def my_annotate(ax, annotation, xy, fontsize=9, rotation=0):
    ax.annotate(annotation, xy=(0,0), xytext=xy, textcoords='axes fraction',
                va='center', ha='center', fontsize=fontsize,
                annotation_clip=False, rotation=rotation, zorder=13)

my_annotate(ax, 'NA', (0.5, -0.10), fontsize=12)
my_annotate(ax, 'Angle Between Arms', (-0.125, 0.5), fontsize=12, rotation=90)
my_annotate(ax, 'Objectives collide\nwith cover slip', (0.7, 0.8))
my_annotate(ax, 'Objectives collide\nwith each other', (0.7, 0.20))
my_annotate(ax, 'Objectives\ncollide with\ncover slip and\n each other', (0.89, 0.5))
my_annotate(ax, 'Feasible', (0.3, 0.5))
my_annotate(ax, 'Epi-detection', (0.5, 0.05))

# Calculate a list of points to sample in region
n_pts = 15
pts = np.mgrid[n/n_pts/2:n:n_pts*1j,0:180:n_pts*1j].reshape((2, n_pts**2)).T.tolist()
        
def is_feasible(pt):
    if pt[1] == 0:
        return True
    elif pt[1] < np.rad2deg(np.pi - 2*np.arcsin(pt[0]/n)) + 8 and pt[1] > np.rad2deg(2*np.arcsin(pt[0]/n)) - 8 and pt[0] < NA_max + 0.1:
        return True
    else:
        return False
    
pts_list = [pt for pt in pts if is_feasible(pt)]
pts = np.array(pts_list).T

# Calculate D for each point
def calc_dispersion(param):
    na = param[0]
    angle = param[1]
    m = multiframe.TwoViewPolScope(n_pts=100, bfp_n=4, 
                                  illum_det_angle=np.deg2rad(angle),
                                  na1=na, na2=na, n_samp=n,
                                  det_type='lens', n_frames=4, 
                                  dist_type='poisson', max_photons=500)
    return util.dispersion_index(m.sa_uncert)

data = []
for i, pt in enumerate(pts.T):
    print('Calculating microscope '+str(i+1)+'/'+str(pts.shape[1]))
    data.append(calc_dispersion(pt))

# Calculate colors
color_map='coolwarm'
color_norm='log'
color_min=1e-3
color_max=1e2
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
sc = ax.scatter(pts[0,:], pts[1,:], c=data, s=0, cmap=cmap, norm=norm, marker='s', lw=0)

# Plot patches
width = n/(n_pts)
for i, (pt, c) in enumerate(zip(pts_list, colors)):
    if pt[1] == 0:
        height = 180/14.5
    if pt[1] == 0:
        height = 180/(n_pts-0.5)
    ax.add_patch(patches.Rectangle((pt[0] - width/2, pt[1] - height/2), width, height, facecolor=c, edgecolor=c))
f.colorbar(sc, cax=cax, orientation='vertical')

# Mask around lines
ax.fill_between(NA, lens_bound, height/2, where=(NA>0.072), color='white', zorder=10)
ax.fill_between(NA, cover_bound, 180, color='white', zorder=10)

print('Saving final figure.')    
f.savefig('sweep_symmetric.pdf', dpi=250)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
