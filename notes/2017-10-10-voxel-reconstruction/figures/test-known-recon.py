from dipsim import multiframe, util, fluorophore, reconstruction
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec
import h5py
import plot_tpc

# Generate truth orientations
nx, ny = 50, 50
t = np.linspace(0, np.pi/2, nx)
p = np.linspace(0, np.pi, ny)
tv, pv = np.meshgrid(t, p)
tvv = np.expand_dims(tv, 2)
pvv = np.expand_dims(pv, 2)
cvv = 3*np.ones((nx, ny, 1))
data = np.concatenate([tvv, pvv, cvv], 2)

# Save truth data
h5f = h5py.File('truth.h5', 'w')
h5f.create_dataset('result', data=data)
h5f.close()

# Plot truth orientations
# plot_tpc.plot_tpc(filename='truth.h5', output='truth.pdf')
# print('X')
# import pdb; pdb.set_trace()

# Create microscope
n = 1.33
alpha1 = 60
na1 = 1.1
na2 = 0.71
theta1 = np.deg2rad(45-12)
theta2 = -np.deg2rad(45+12)
dose_ratio = 3
total_photons = 10000
dose2 = total_photons/(1 + dose_ratio)
dose1 = total_photons - dose2

exp = multiframe.MultiFrameMicroscope(ill_thetas=[theta1, theta2], det_thetas=[theta2, theta1],
                                      ill_nas=[na1, na2], det_nas=[na2, na1],
                                      ill_types=2*['sheet'], det_types=2*['lens'],
                                      colors=['(1,0,0)', '(0,0,1)'], n_frames=4,
                                      max_photons=[dose1, dose2], n_samp=1.33)

# Generate the intensities
def intensities_wrapper(data, exp=exp):
    flu = fluorophore.Fluorophore(*data)
    return exp.calc_intensities(fluorophore=flu)

intensities = np.apply_along_axis(intensities_wrapper, axis=2, arr=data)
import pdb; pdb.set_trace()

# Perform reconstructions
def recon_wrapper(data, multiframe=exp):
    recon = reconstruction.Reconstruction(multiframe=multiframe, recon_type='tpc',
                                          data=data)
    recon.evaluate()
    result = recon.estimated_fluorophore.theta, recon.estimated_fluorophore.phi, recon.estimated_fluorophore.c
    print(data, result)    
    return result

result = np.apply_along_axis(recon_wrapper, axis=2, arr=intensities, multiframe=exp)
import pdb; pdb.set_trace()

# Save result
h5f = h5py.File('data.h5', 'w')
h5f.create_dataset('result', data=result)
h5f.close()

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
