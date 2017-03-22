from dipsim import fluorophore, illuminator, detector, multiframe, microscope, stats, util
import numpy as np
import matplotlib.pyplot as plt
import time; start = time.time(); print('Running...')

m = multiframe.NFramePolScope(n_frames=1)
m.plot_solid_angle_min_std('random.png', 'Random', n=2000, display='interact', random_colors=True, show_edges=True)

print('Total time:', np.round(time.time() - start, 2), 's')
