from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time;  print('Running...')

# Compute and plot on axes
pol = np.pi/4
det_axis = np.pi/4
ill_axis = -np.pi/4
det_na = 0.8
ill_na = 0.8
kappa = 1

# Create microscope
ill = illuminator.Illuminator(illum_type='wide',
                              theta_optical_axis=ill_axis,
                              na=ill_na,
                              phi_pol=pol)

det = detector.Detector(theta_optical_axis=np.array(det_axis),
                        na=det_na)

m = microscope.Microscope(illuminator=ill, detector=det, max_photons=100)

# Calculate total efficiency
epsrel = 1000
kappas = [None, 0]#np.arange(-10, 10, 1)
for kappa in kappas:
    start = time.time();
    x = m.calc_sensitivity(direction=(0,0), kappa=kappa, epsrel=epsrel)
    print('kappa = '+str(kappa)+'\t'+str(np.round(x, 5))+'\t'+str(np.round(time.time() - start, 5)))
#os.system('say "done"')
