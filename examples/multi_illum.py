from dipsim import fluorophore, illuminator, detector, multiframe, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

for i in np.arange(2, 3):
    m = multiframe.NFramePolScope(n_frames=i)
    m.plot_solid_angle_min_std('solid'+str(i)+'.png', n=5000, display='interact')

print('Total time:', np.round(time.time() - start, 2), 's')
