from dipsim import fluorophore, illuminator, detector, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

for bfp_rad in np.arange(1, 7, 1):
    f = 10
    bfp_n = 300

    illx = illuminator.Illuminator(illum_type='kohler',
                                   optical_axis=np.array([0., 0., 1.]),
                                   f=f,
                                   bfp_rad=bfp_rad,
                                   bfp_pol=np.array([1., 0., 0.]),
                                   bfp_n=bfp_n)

    det = detector.Detector(optical_axis=np.array([0, 0, 1]),
                            na=1.5,
                            n=1.5)

    m = microscope.Microscope(illuminator=illx, detector=det)

    m.plot_intensities_from_single_fluorophore('intensity'+str(bfp_rad)+'.png', 'Intensity', n=5000, display='save')

