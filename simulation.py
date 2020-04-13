#!/usr/bin/env python3.6

import numpy as np

class trace_generator:
    '''generating traces of two levels drawn from poisson distribution with one changepoint'''
    def __init__(self,
                time_constant   = 400,
                time_end        = 800,
                n_points        = 50,
                n_traces        = 2000,
                bright_level    = 20,
                dark_level      = 8):

        self.time_constant  = time_constant
        self.time_end       = time_end
        self.n_points       = n_points
        self.n_traces       = n_traces
        self.bright_level   = bright_level
        self.dark_level     = dark_level

    def draw_loss_times(self):
        '''draw times from exponential distribution'''

        size    = self.n_traces
        times   = []
        while len(times) < size:
            #print(len(times))

            time = np.random.exponential(scale = self.time_constant)

            if time < 0.9*self.time_end and time > 0.1*self.time_end:
                times.append(time)

        return times

    def generate_traces(self,times):
        '''draw from poisson'''

        traces = []
        for time in times:

            bright_length   = int(self.n_points*time/self.time_end)
            bright_traces   = np.random.poisson(lam = self.bright_level, size = bright_length)
            dark_length     = self.n_points - bright_length

            if dark_length > 0:
                dark_traces     = np.random.poisson(lam = self.dark_level, size = dark_length)
                traces.append(np.concatenate((bright_traces,dark_traces)))

        return traces
