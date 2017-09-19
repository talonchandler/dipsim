from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import scipy.special as special


inch_fig = 3
f, axs = plt.subplots(nrows=1, ncols=7, figsize=(7*inch_fig, inch_fig), subplot_kw={'projection':'3d'})
plt.subplots_adjust(wspace=-0.7)
kappas = [-3, -2, -.5, 0, .5, 2, 3]
for ax, kappa in zip(axs, kappas):
    ax.set_aspect("equal")
    
    phi, theta = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
    norm = 1.0/(4*np.pi*special.hyp1f1(0.5, 1.5, kappa))
    weight = np.exp(kappa*(np.cos(theta)**2))
    x = norm*weight*np.cos(phi)*np.sin(theta)
    y = norm*weight*np.sin(phi)*np.sin(theta)
    z = norm*weight*np.cos(theta)
    ax.plot_wireframe(x, y, z, color="k", lw=0.1)

    ax.set_xlim(-0.3, 0.3)
    ax.set_ylim(-0.3, 0.3)
    ax.set_zlim(-0.3, 0.3)    
    ax.set_axis_off()
    ax.azim = 0
    ax.elev = 30
    ax.set_title('$\kappa = '+str(kappa)+'$', y=0.9, fontsize=20)
f.savefig('watson.pdf')
