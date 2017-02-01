import time; start = time.time(); print('Running...')
from dipsim import dipole, exciter, detector, microscope
import numpy as np
import matplotlib.pyplot as plt

dip = dipole.Dipole(position=np.array([0, 0, 0]),
                    orientation=np.array([1, 1, 1]),
                    photon_yield=1)

e = exciter.Exciter(optic_axis=np.array([0, 0, 0]))

det = detector.Detector(k=2*np.pi/(600e-9), na=1.0, n0=1.518, n1=1.0,
                        optic_axis=np.array([0, 0, 1]),
                        img_axis=np.array([1, 0, 0]),
                        f_obj=3e-3,
                        n_pixel=2**8,
                        d_pixel=2./(2**8))

# det = detector.Detector(k=1, n1=1,
#                         optic_axis=np.array([1, 0, 0]),
#                         img_axis=np.array([0, 0, 1]),
#                         f_obj=1,
#                         n_pixel=256,
#                         d_pixel=10./256)

m = microscope.Microscope(exciter=e,
                          detector=det)

start2 = time.time(); 
m.calc_bfp(n_phi=256, n_theta=64, dipole=dip)
m.calc_img()
print('Calc time:', np.round(time.time() - start2, 2), 's')

m.plot_detector('z_det1.png')

m.detector.optic_axis=np.array([1., 0, 0])
m.detector.img_axis=np.array([0, 0, 1.])
m.calc_bfp(n_phi=256, n_theta=64, dipole=dip)
m.plot_detector('x_det1.png')

m.detector.optic_axis=np.array([1., 1., 1.])/np.sqrt(3)
m.detector.img_axis=np.array([1.,-1.,0])/np.sqrt(2)
m.calc_bfp(n_phi=256, n_theta=64, dipole=dip)
m.plot_detector('oblique_det1.png')


m.detector.optic_axis=np.array([0, 0, 1.])
m.detector.img_axis=np.array([1., 0, 0])
dip.orientation = np.array([0, 0, 1])
m.calc_bfp(n_phi=256, n_theta=64, dipole=dip)
m.plot_detector('z_det2.png')

m.detector.optic_axis=np.array([1., 0, 0])
m.detector.img_axis=np.array([0, 0, 1.])
m.calc_bfp(n_phi=256, n_theta=64, dipole=dip)
m.plot_detector('x_det2.png')


print('Total time:', np.round(time.time() - start, 2), 's')

