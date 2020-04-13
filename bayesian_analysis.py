#!/usr/bin/env python3.6

import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import poisson
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.misc import factorial
from scipy.special import gamma


class bayesian:
    def __init__(self):
        pass

    def prob_trace(self,counts,r1,r2):
        '''takes in one trace and outputs prob trace'''
        prob_lost   = []
        T           = len(counts)
        for i in range(T):

            N_photons_before            = np.sum(counts[0:i+1])
            N_photons_after             = np.sum(counts[i:])
            expected_N_photons_before   = r1*i
            expected_N_photons_after    = r2*(T-i)

            p1      = poisson.pmf(N_photons_before, expected_N_photons_before)
            p2      = poisson.pmf(N_photons_after,expected_N_photons_after)
            prob    = p1*p2

            prob_lost.append(prob)

        return prob_lost

    def double_poisson(self,k, *params):
        (r1,r2,p1,p2)   = params
        poisson_1       = p1 * (r1**k/factorial(k)) * np.exp(-r1)
        poisson_2       = p2 * (r2**k/factorial(k)) * np.exp(-r2)
        return  poisson_1 + poisson_2

    def single_poisson(self,k, lam):
        poisson         = (lam**k/factorial(k)) * np.exp(-lam)
        return  poisson

    def fit_hist(self,traces, display = False):
        '''takes in all traces and returns the means of the double histogram'''
        bins_for_fit    =   40

        all_counts      = np.ndarray.flatten(np.array(traces))
        points, edges   = np.histogram(all_counts,bins_for_fit)

        x_axis          = [(edges[i]+edges[i+1])/2.0 for i in range(len(edges)-1)]
        all_mean        = np.mean(all_counts)
        p0              = [all_mean*.5,all_mean*1.5,0.8,0.2]

        popt, pcov      = curve_fit(self.double_poisson,x_axis,
                                points,p0=p0)

        bright_cnt      = max(popt[0],popt[1])
        dark_cnt        = min(popt[0],popt[1])
        percent_bright  = popt[2]/(popt[2]+popt[3])
        percent_dark    = 1-percent_bright

        print("avg bright/dark counts = ", bright_cnt, dark_cnt)
        print("percent bright/dark = ", percent_bright, percent_dark)

        if display:
            plt.figure()
            plt.hist(all_counts, bins = bins_for_fit,
                    color="blue", alpha = 0.5)

            x_plot  = np.linspace(0, max(x_axis), bins_for_fit*100)
            hist1   = popt[2]*self.single_poisson(x_plot, popt[0])
            hist2   = popt[3]*self.single_poisson(x_plot, popt[1])

            plt.plot(x_plot, hist1, 'b.')
            plt.plot(x_plot, hist2, 'b.-')
            plt.xlabel("counts per 16 us")
            plt.ylabel("frequency")
            plt.show()

        return bright_cnt, dark_cnt


    def prior(lam,k,theta):
        return lam**(k-1)*np.exp(lam/theta)

    def log_g_over_s(self,r,s):
        e = np.exp(1)

        if r < 50:
            return np.log(gamma(r)/s**r)
        else:
            ans  = r*np.log(r/e/s)
            ans += 0.5*np.log(2*np.pi/r)
            return ans

    #def bayes_factor(n,)


    def return_normed_prob(self,log_prob):
        max_log     = max(log_prob)
        log_norm    = np.subtract(log_prob, max_log)
        lin_prob    = []

        for logp in log_norm:
            if logp < -150:
                lin_prob.append(0)
            else:
                lin_prob.append(np.exp(logp))

        return lin_prob


    def changepoint(self,counts):
        '''takes in one trace and outputs prob trace'''

        #prior settings
        k = 0.5; theta = 1e6;

        log_probs       = []
        T               = len(counts)
        bay_sum         = 0

        for i in range(T):
            N_before    = np.sum(counts[0:i+1])
            N_after     = np.sum(counts[i:])
            r1          = N_before + k
            s1          = i+1e-4
            r2          = N_after + k
            s2          = T-i+1e-4

            log_p       = self.log_g_over_s(r1,s1)
            log_p      += self.log_g_over_s(r2,s2)

            log_probs.append(log_p)

            #bay_sum    += np.exp(log_p-np.log(gamma(T+0.5)))

        log_probs       = log_probs[1:]
        lin_probs       = [0]+self.return_normed_prob(log_probs)

        #bayes_factor    = 4*np.sqrt(np.pi)/bay_sum

        #print("bayes factor", bayes_factor)

        return np.multiply(lin_probs,1/T)
