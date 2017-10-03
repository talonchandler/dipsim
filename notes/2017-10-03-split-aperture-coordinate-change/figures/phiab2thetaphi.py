import numpy as np
import util
import matplotlib.pyplot as plt
import os

def fibonacci_sphere(n):
    # Returns "equally" spaced points on a unit sphere in spherical coordinates.
    # http://stackoverflow.com/a/26127012/5854689
    z = np.linspace(1 - 1/n, -1 + 1/n, num=n) 
    theta = np.arccos(z)
    phi = np.mod((np.pi*(3.0 - np.sqrt(5.0)))*np.arange(n), 2*np.pi) - np.pi
    return np.vstack((theta, phi)).T

def phi_prime(theta, phi, alpha=np.pi/6):
    # Returns phi coordinate in a frame rotated by a right-handed rotation by
    # angle alpha about the positive y-axis. 
    if theta == 0:
        return 0
    num = np.cos(alpha)*np.cos(phi)*np.sin(theta) - np.sin(alpha)*np.cos(theta)
    den = np.sqrt(1 - (np.sin(alpha)*np.cos(phi)*np.sin(theta) + np.cos(alpha)*np.cos(theta))**2)
    if phi < np.pi and phi > 0:
        return np.arccos(num/den)
    elif phi < 0 and phi > -np.pi:
        return -np.arccos(num/den)
    else:
        return 0

def thetaphi2phiab(theta, phi, alpha=np.pi/6):
    # Converts theta, phi coordinates into phi_a, phi_b coordinates
    phi_a = phi_prime(theta, phi, alpha=alpha)
    phi_b = phi_prime(theta, phi, alpha=-alpha)
    return np.array([phi_a, phi_b])

def generate_LUTs(n=1e4, alpha=np.pi/6):
    # Generates LUTs in both coordinate systems
    thetaphiLUT = fibonacci_sphere(n)
    phiabLUT = []
    for direction in thetaphiLUT:
        phiabLUT.append(thetaphi2phiab(*direction, alpha=alpha))
    phiabLUT = np.array(phiabLUT)
    return thetaphiLUT, phiabLUT

def phiab2thetaphi(phia, phib, thetaphiLUT=None, phiabLUT=None):
    # Converts phi_a, phi_b coordinates into theta, phi using lookup tables.
    if phiabLUT is None:
        thetaphiLUT, phiabLUT = generate_LUTs()
    target = np.array([phia, phib])
    min_ind = np.argmin(np.linalg.norm(target - phiabLUT, axis=1))
    return thetaphiLUT[min_ind]

# Test error introduced by reconstruction
thetaphiLUT, phiabLUT = generate_LUTs(n=5e5)
import os; import time; start = time.time(); print('Running...')
print(phiab2thetaphi(0, 0.1, thetaphiLUT=thetaphiLUT, phiabLUT=phiabLUT))
print('Total time: '+str(np.round(time.time() - start, 4)))

directions = fibonacci_sphere(20000)
results = []
for direction in directions:
    noisy_phiab = thetaphi2phiab(*direction) + np.random.normal(0, np.deg2rad(0))
    recon_direction = phiab2thetaphi(*noisy_phiab, thetaphiLUT=thetaphiLUT, phiabLUT=phiabLUT)
    results.append(np.abs(recon_direction - direction))

results = np.array(results)

# Plots results
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
plt.subplots_adjust(wspace=0.1, hspace=-0.1)
caxs = util.generate_caxs(axs)

util.plot_sphere('', directions=directions, data=np.rad2deg(results[:,0]), my_ax=axs[0], my_cax=caxs[0], color_norm='log',
                 color_min=1e-2, color_max=1e2)
util.plot_sphere('', directions=directions, data=np.rad2deg(results[:,1]), my_ax=axs[1], my_cax=caxs[1], color_norm='log',
                 color_min=1e-2, color_max=1e2)

axs[0].annotate('$\\theta$ error [degrees]', xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction',
                va='center', ha='center', rotation=0, fontsize=16, annotation_clip=False)
axs[1].annotate('$\\phi$ error [degrees]', xy=(0,0), xytext=(0.5, 1.1), textcoords='axes fraction',
                va='center', ha='center', rotation=0, fontsize=16, annotation_clip=False)

fig.savefig('recon-error.pdf', dpi=300)
os.system('say "done"')
