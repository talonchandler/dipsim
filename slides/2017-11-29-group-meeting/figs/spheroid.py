from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import scipy.special as special


inch_fig = 3
f, axs = plt.subplots(nrows=1, ncols=7, figsize=(7*inch_fig, inch_fig), subplot_kw={'projection':'3d'})
plt.subplots_adjust(wspace=-0.7)
kappas = [0.1, 0.25, 0.5, 1, 2, 4, 10]
for ax, kappa in zip(axs, kappas):
    ax.set_aspect("equal")
    
    phi, theta = np.mgrid[0:2*np.pi:15j, 0:np.pi:20j]
    if kappa == 1:
        norm = 4*np.pi
    elif kappa < 1:
        z = np.sqrt(1 - kappa**2)
        norm = 2*np.pi*(1 + (kappa**2/(2*z))*np.log((1+z)/(1-z)))
    elif kappa > 1:
        z = np.sqrt(kappa**2 - 1)        
        norm = 2*np.pi*(1 + (kappa**2/z)*np.arcsin(z/kappa))
    weight = 1.5
    x = (1.0/norm)*weight*np.cos(phi)*np.sin(theta)
    y = (1.0/norm)*weight*np.sin(phi)*np.sin(theta)
    z = (1.0/norm)*weight*np.cos(theta)*kappa
    ax.plot_wireframe(x, y, z, color="k", lw=0.1)

    ax.set_xlim(-0.3, 0.3)
    ax.set_ylim(-0.3, 0.3)
    ax.set_zlim(-0.3, 0.3)    
    ax.set_axis_off()
    ax.azim = 0
    ax.elev = 30
    ax.set_title('$\kappa = '+str(kappa)+'$', y=0.9, fontsize=20)
f.savefig('spheroid.pdf')
