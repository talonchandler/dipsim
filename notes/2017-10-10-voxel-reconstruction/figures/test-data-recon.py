from dipsim import multiframe, util, fluorophore, reconstruction
import numpy as np
import matplotlib.pyplot as plt
import os; import time; start = time.time(); print('Running...')
import matplotlib.gridspec as gridspec
import h5py

# Import data
name_head = '/Users/Talon/Dropbox/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['B_reg_P3.tif', 'B_reg_P4.tif', 'B_reg_P1.tif', 'B_reg_P2.tif',
         'A_reg_P3.tif', 'A_reg_P4.tif', 'A_reg_P1.tif', 'A_reg_P2.tif']
input_files = [name_head + name for name in names]
data = np.zeros((50, 50, 8))
for i, input_file in enumerate(input_files):
    im = util.tiff2array(input_file, x=215, y=585, z=175, width=50, height=50, slices=1)
    data[:,:,i] = im

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

# Perform reconstructions
def recon_wrapper(data, multiframe=exp):
    recon = reconstruction.Reconstruction(multiframe=multiframe, recon_type='tpc',
                                          data=data)
    recon.evaluate()
    theta_prime = util.theta_prime(recon.estimated_fluorophore.theta, recon.estimated_fluorophore.phi, -np.deg2rad(45-12))
    phi_prime = util.phi_prime(recon.estimated_fluorophore.theta, recon.estimated_fluorophore.phi, -np.deg2rad(45-12))    
    result = theta_prime, phi_prime, recon.estimated_fluorophore.c
    
    print(data, result)    
    return result

result = np.apply_along_axis(recon_wrapper, axis=2, arr=data, multiframe=exp)

# Save result
h5f = h5py.File('data.h5', 'w')
h5f.create_dataset('result', data=result)
h5f.close()

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
