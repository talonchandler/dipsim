from dipsim import fluorophore, illuminator, detector, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

for bfp_rad in np.arange(4, 7, 1):
    f = 10
    bfp_n = 100
    illy = illuminator.Illuminator(illum_type='kohler',
                                   optical_axis=np.array([0., 0., 1.]),
                                   f=f,
                                   bfp_rad=bfp_rad,
                                   bfp_pol=np.array([0., 1., 0.]),
                                   bfp_n=bfp_n)

    illx = illuminator.Illuminator(illum_type='kohler',
                                   optical_axis=np.array([0., 0., 1.]),
                                   f=f,
                                   bfp_rad=bfp_rad,
                                   bfp_pol=np.array([1., 0., 0.]),
                                   bfp_n=bfp_n)

    illxy = illuminator.Illuminator(illum_type='kohler',
                                   optical_axis=np.array([0., 0., 1.]),
                                   f=f,
                                   bfp_rad=bfp_rad,
                                   bfp_pol=np.array([1./np.sqrt(2), 1./np.sqrt(2), 0.]),
                                   bfp_n=bfp_n)
    
    print(bfp_rad, illx.E_eff, illy.E_eff, illxy.E_eff)
