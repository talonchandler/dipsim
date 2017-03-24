from dipsim import fluorophore, illuminator, detector, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

n = 100
dpi = 50
vis_px = 750

f = 10
bfp_n = 16
bfp_rad = 1
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

m.plot_intensities_from_single_fluorophore('intensity'+str(bfp_rad)+'.png',
                                           n=n, interact=False, color_norm='linear',
                                           dpi=dpi, vis_px=vis_px)

m.plot_intensities_from_single_fluorophore_new('intensity_new'+str(bfp_rad)+'.png',
                                               n=n, interact=False, color_norm='linear',
                                               dpi=dpi, vis_px=vis_px)

m.draw_scene('scene'+str(bfp_rad)+'.png', interact=False,
             dpi=dpi, vis_px=vis_px)
print('Total time: '+str(np.round(time.time() - start, 2)))
