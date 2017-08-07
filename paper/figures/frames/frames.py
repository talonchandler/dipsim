import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
import subprocess    

inch_fig = 3.3
fig_labels = ['a)', 'b)', 'c)']

# Setup figure and axes
fig = plt.figure(figsize=(3*inch_fig, 1*inch_fig))
gs0 = gridspec.GridSpec(1, 3, wspace=0, hspace=0)
ax0 = plt.subplot(gs0[0])
ax1 = plt.subplot(gs0[1])
ax2 = plt.subplot(gs0[2])

for ax, fig_label  in zip([ax0, ax1, ax2], fig_labels):
    ax.annotate(fig_label, xy=(0,0), xytext=(0, 0.95), textcoords='axes fraction',
                va='bottom', ha='center', fontsize=14, annotation_clip=False)

def plot_img(name, ax):
    subprocess.call(['asy', name+'.asy'])
    subprocess.call(['convert', '-density', str(300), '-units', 'PixelsPerInch', name+'.pdf', name+'.png'])
    im = mpimg.imread(name+'.png')
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticklabels([])
    ax.imshow(im, interpolation='none')


plot_img('frame_a', ax0)
plot_img('frame_b', ax1)
plot_img('frame_c', ax2)
    
print('Saving final figure.')    
fig.savefig('../../paper/frames.pdf', dpi=250)
    
