import matplotlib.pyplot as plt
import numpy as np
from dipsim import util

# Choose files
name_head = '/Users/Talon/Dropbox/20170725_Bob_Actin_results/Cell1_LSimaging_registerred/SPIM'
names = ['A_reg_P1.tif', 'A_reg_P2.tif', 'A_reg_P3.tif', 'A_reg_P4.tif',
         'B_reg_P1.tif', 'B_reg_P2.tif', 'B_reg_P3.tif', 'B_reg_P4.tif']
input_files = [name_head + name for name in names]
input_files = np.array(input_files).reshape(2, 4)

# Choose labels
col_labels = [0, 45, 90, 135]
col_labels = ['$'+str(lab)+'^{\circ}$' for lab in col_labels]
row_labels = ['View A: 1.1 NA', 'View B: 0.71 NA']

# Plot array
util.plot_array(input_files, 'octet.pdf',
                row_labels=row_labels, col_labels=col_labels,
                line_start=(215, 635), line_end=(265, 585), zslice=175, 
                roi_upper_left=(215, 585), roi_wh=(50, 50))

util.plot_array(input_files, 'roi-octet.pdf',
                row_labels=row_labels, col_labels=col_labels,
                line_start=(0, 49), line_end=(49, 0), zslice=175,
                roi_upper_left=(215, 585), roi_wh=(50, 50), plot_roi=True)
