#!/usr/bin/env python3.6

import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import poisson
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import factorial

def adjust_spines(ax,spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
        else:
            spine.set_color('none') # don't draw spine

    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        ax.xaxis.set_ticks([])


class edge_detection:
    def __init__(self):
        pass

    def smoothen_and_convolve(self, time, trace, display = False):
        '''takes in one trace. calculates the slope of its moving
        average. N is chosen for good signal to noise, but we also
        suffer from the edge effect'''

        N   = 5#int(0.05*len(trace))
        smooth_trace = np.convolve(trace, np.ones((N,))/N, mode='same')
        slope_trace  = np.convolve(smooth_trace,[1,0,0,-1], mode='same')

        return smooth_trace, slope_trace

    def plot_three_traces(self, time, trace, smooth_trace, slope_trace):
        fig = plt.figure(figsize=(10,6))
        ax  = fig.add_subplot(3,1,1)

        ax.plot(time,trace,color="r",marker='o',linestyle='--')
        ax.set_ylabel("counts per bin")
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        adjust_spines(ax,['left'])

        ax  = fig.add_subplot(3,1,2)

        ax.plot(time,smooth_trace,color = "#003F87",marker='o',linestyle='-.')
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))

        ax.set_ylabel("moving average (a.u.)")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        adjust_spines(ax,['left'])

        ax          = fig.add_subplot(3,1,3)
        ax.plot(time[5:-5],slope_trace[5:-5],color = "#003F87",marker='o',linestyle='-.')
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))

        ax.set_ylabel("slope (a.u.)")
        ax.set_xlabel("time [us]")
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        adjust_spines(ax,['left','bottom'])

        plt.show()
