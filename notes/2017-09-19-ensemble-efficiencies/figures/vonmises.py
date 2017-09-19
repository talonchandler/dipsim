from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import scipy.special as special

inch_fig = 3
f, axs = plt.subplots(nrows=1, ncols=7, figsize=(7*inch_fig, inch_fig), subplot_kw={'projection':'3d'})
plt.subplots_adjust(wspace=-0.7)
kappas = [0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0]
for ax, kappa in zip(axs, kappas):
    ax.set_aspect("equal")
    
    phi, theta = np.mgrid[0:2*np.pi:20j, 0:np.pi:30j]
    norm = np.sqrt(kappa)/(((2*np.pi)**(1.5))*special.iv(0.5, kappa))
    weight = np.exp(kappa*np.cos(theta))
    x = norm*weight*np.cos(phi)*np.sin(theta)
    y = norm*weight*np.sin(phi)*np.sin(theta)
    z = norm*weight*np.cos(theta)
    ax.plot_wireframe(x, y, z - 0.7, color="k", lw=0.1)

    ax.set_xlim(-0.9, 0.9)
    ax.set_ylim(-0.9, 0.9)
    ax.set_zlim(-0.9, 0.9)    
    ax.set_axis_off()
    ax.azim = 0
    ax.elev = 30
    ax.set_title('$\kappa = '+str(kappa)+'$', y=0.9, fontsize=20)
f.savefig('vonmises.pdf')
