from dipsim import fluorophore, illuminator, detector, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')
import functools

ill = illuminator.Illuminator(illum_type='kohler',
                              optical_axis=np.array([0., 0., 1.]),
                              f=1.0,
                              bfp_rad=0.1,
                              bfp_pol=np.array([0., 1., 0.]),
                              bfp_n=64)

det = detector.Detector(optical_axis=np.array([0, 0, 1]),
                        na=1.5,
                        n=1.5)

m = microscope.Microscope(illuminator=ill,
                          detector=det)

def my_micro(m, args):
    theta = args[0]
    phi = args[1]
    
    flu_dir = np.array([np.sin(theta)*np.cos(phi),
                       np.sin(theta)*np.sin(phi),
                       np.cos(theta)])
    
    flu = fluorophore.Fluorophore(position=np.array([0, 0, 0]),
                                  mu_abs=flu_dir,
                                  mu_em=flu_dir)

    I = m.calc_total_intensity([flu])
    
    n_photons = 1e6
    return n_photons*I

my_pdf = functools.partial(my_micro, m)

pdf = stats.Pdf(my_pdf, dist_type='poisson')
#crlb = pdf.crlb([0.1, 2], [0.00001, 0.00001])

# Calculate std_solid_angle
def std_solid_angle(x0, dx):
    var = pdf.crlb(x0, dx)
    return np.sin(x0[0])*np.sqrt(var[0])*np.sqrt(var[1])

# Calculate for each point
n = 100
tt, pp = np.meshgrid(np.linspace(0, np.pi, n), np.linspace(0, 2*np.pi, n))
directions = np.stack((tt, pp))
crlb_all = np.apply_along_axis(std_solid_angle, 0, directions, [0.00001, 0.00001])
crlb = np.apply_along_axis(pdf.crlb, 0, directions, [0.00001, 0.00001])
intensity = np.apply_along_axis(my_pdf, 0, directions) 

util.plot_sphere('solid.png', tt, pp, crlb_all)
util.plot_sphere('theta.png', tt, pp, np.sqrt(crlb[0, :, :]))
util.plot_sphere('phi.png', tt, pp, np.sqrt(crlb[1, :, :]))
util.plot_sphere('intensity.png', tt, pp, intensity)

print('Total time:', np.round(time.time() - start, 2), 's')


