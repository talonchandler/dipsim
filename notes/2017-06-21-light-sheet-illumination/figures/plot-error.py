from dipsim import util
import numpy as np
import matplotlib.pyplot as plt
from functools import partial

import matplotlib
matplotlib.rcParams.update({'font.size': 14})

# All lengths in um

# Define functions
def w(z, z0, w0):
    return w0*np.sqrt(1+(z/z0)**2)

def z_excite(x, z, w0, phi_pol, lamb=0.488):
    k = 2*np.pi/lamb
    z0 = k*(w0**2)/2.0
    a = w0/w(z, z0, w0)
    b = np.exp(-(2*x**2)/(w(z, z0, w0)**2))
    c = ((x*np.cos(phi_pol))**2 + (0.5*w(z, z0, w0)*np.sin(phi_pol))**2)/(z0**2)
    return a*b*c

# Setup x-z plots
n_rows = 1
n_cols = 2
inch_fig = 5
size = (inch_fig*n_cols, inch_fig*n_rows)
fig, axs = plt.subplots(n_rows, n_cols, figsize=size)
plt.subplots_adjust(wspace=0.5, hspace=0)
caxs = util.generate_caxs(axs)

x = np.arange(-40, 40, 0.05)
z = np.arange(-40, 40, 0.05)
X, Z = np.meshgrid(x, z)

phi_pols = [0, np.pi/2]
for i, phi_pol in enumerate(phi_pols):
    I = z_excite(X, Z, w0=1.2, phi_pol=phi_pol, lamb=0.488)
    extent = (x[0], x[-1], z[0], z[-1])
    cmap = axs[i].imshow(I, interpolation='none', cmap='Greys', extent=extent)
    fig.colorbar(cmap, cax=caxs[i], orientation='vertical')

    axs[i].set_xlabel('$x$ ($\mu$m)')
    axs[i].set_ylabel('$z$ ($\mu$m)')
    if i == 0:
        axs[i].set_title('$\phi_{\mathrm{pol}} = 0$', y=1.05)
    if i == 1:
        axs[i].set_title('$\phi_{\mathrm{pol}} = \pi/2$', y=1.05)

fig.savefig('error.pdf')

# Setup max error plots
def find_max_error(waist, FOV):
    x = np.arange(-(FOV/2.0), (FOV/2.0), 0.1)
    z = np.arange(-(FOV/2.0), (FOV/2.0), 0.1)
    X, Z = np.meshgrid(x, z)
    I = z_excite(X, Z, w0=waist, phi_pol=np.pi/2, lamb=0.488)
    return np.max(I)

fig, ax = plt.subplots(1, 1, figsize=(5, 5))

FOVs = [80]
for FOV in FOVs:
    find_max_error_partial = partial(find_max_error, FOV=FOV)
    find_max_error_vector = np.vectorize(find_max_error_partial)
    w0 = np.arange(1, 5.1, 0.1)
    max_error = find_max_error_vector(w0)
    ax.plot(w0, max_error, '-k')

ax.set_yscale('log')
ax.set_xlim([1,5])
ax.set_xlabel('$w_0$ ($\mu$m)')
ax.set_ylabel('Maximum Intensity Weighted Excitation Ratio')
fig.savefig('max-error.pdf')
