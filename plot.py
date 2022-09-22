import matplotlib.pyplot as plt
import numpy as np

def set_plot(ratio=(13,9), font_size=22, ticks=[0,1], ticks_step=0.1, labels=['x','y', 'title']):
    plt.figure(figsize=ratio)
    plt.rcParams.update({'font.size': font_size})
    plt.xlabel(labels[0], fontsize=font_size)
    plt.ylabel(labels[1], fontsize=font_size)
    plt.title(labels[2], fontsize=font_size)
    #--------ticks
    plt.xticks(np.arange(ticks[0], ticks[1]+1, step=ticks_step))
    plt.tick_params(axis='both', which='major', labelsize=font_size)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #--------plot
    plt.grid()