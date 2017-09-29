from dipsim import multiframe, util, detector, illuminator, microscope, util
import numpy as np
import matplotlib.pyplot as plt

import os; import time; start = time.time(); print('Running...')

print('Computing microscope: ')
# y-polarized input from x optical axis
ill = illuminator.Illuminator(illum_type='sheet', theta_optical_axis=np.pi/2, phi_pol=np.pi/2)
# perp: x-polarized detection
# par: y-polarized detection
det_perp = detector.Detector(det_type='polarized', na=0.001, phi_pol=0)
det_par = detector.Detector(det_type='polarized', na=0.001, phi_pol=np.pi/2)
m_perp = microscope.Microscope(illuminator=ill, detector=det_perp, max_photons=1e10)
m_par = microscope.Microscope(illuminator=ill, detector=det_par, max_photons=1e10)

f, axs = plt.subplots(1, 2, figsize=(10, 5))
util.draw_scene(m_perp.scene_string(), my_ax=axs[0], dpi=300)
util.draw_scene(m_par.scene_string(), my_ax=axs[1], dpi=300)
axs[0].annotate('$I_\perp$', xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)
axs[1].annotate('$I_\parallel$', xy=(0,0), xytext=(0.5, 1.05), textcoords='axes fraction',
            va='bottom', ha='center', fontsize=14, annotation_clip=False)

f.savefig('schematic.pdf', dpi=300)

print('Total time: '+str(np.round(time.time() - start, 2)))
os.system('say "done"')
