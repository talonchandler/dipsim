from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

# Compute anisotropy
def anisotropy(phi):
    print('Computing microscope: ')
    # y-polarized input from x optical axis
    ill = illuminator.Illuminator(illum_type='sheet', theta_optical_axis=np.pi/2, phi_pol=np.pi/2)
    # perp: x-polarized detection
    # par: y-polarized detection
    det_perp = detector.Detector(det_type='polarized', na=0.001, phi_pol=np.pi/2)
    det_par = detector.Detector(det_type='polarized', na=0.001, phi_pol=0)
    m_perp = microscope.Microscope(illuminator=ill, detector=det_perp, max_photons=1e10)
    m_par = microscope.Microscope(illuminator=ill, detector=det_par, max_photons=1e10)
    I_perp = m_perp.calc_intensity((np.pi/2, phi), kappa=0.0, epsrel=1e-8)
    I_par = m_par.calc_intensity((np.pi/2, phi), kappa=0.0, epsrel=1e-8)
    r = (I_par - I_perp)/(I_par + 2*I_perp)
    print(I_perp, I_par, r)
    if I_par + 2*I_perp == 0:
        print("*")
        return 0.0
    else:
        return r

f, ax = plt.subplots(1, 1, figsize=(5, 5))
x = list(np.linspace(-np.pi, np.pi, 100))
v_aniso = np.vectorize(anisotropy)
r_vect = v_aniso(x) 
ax.plot(x, r_vect, '-k')
ax.set_xlabel('$\phi$')
ax.set_ylabel('$r$')
f.savefig('anisotropy.pdf', dpi=300)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
