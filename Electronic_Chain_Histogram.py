import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from matplotlib.offsetbox import AnchoredText
import matplotlib.mlab as mlab
from scipy.stats import norm
from scipy import stats

from charge_RC import charge


#The purpose of this code is to display and save the histograms plotting the maximum charge extracted from each waveform.

#for proper error handling for bins with counts less than 20, symmetric error bars are pulled from this method
def error_bins(raw_data, fitted_data, raw_error):
    chi2_table = {0: [0.01, 1.29], 1: [0.27, 2.75], 2: [0.74, 4.25], 3: [1.10, 5.30], 4: [2.34, 6.78], 5: [2.75, 7.81],
                  6: [3.82, 9.28],
                  7: [4.25, 10.30], 8: [5.30, 11.32], 9: [6.33, 12.79], 10: [6.78, 13.81], 11: [7.81, 14.82],
                  12: [8.83, 16.29], 13: [9.28, 17.30],
                  14: [10.30, 18.32], 15: [11.32, 19.32], 16: [12.33, 20.80], 17: [12.79, 21.81], 18: [13.81, 22.82],
                  19: [14.82, 23.82], 20: [15.83, 25.30]}

    for i in range(len(raw_data)):
        if raw_data[i] <= 10:
            error_band = chi2_table.get(raw_data[i])
            upper = abs(error_band[0] - raw_data[i])
            lower = abs(error_band[1] - raw_data[i])
            residual = fitted_data[i] - raw_data[i]
            if(residual >= 0):
                raw_error[i] = upper
            elif(residual < 0):
                raw_error[i] = lower
    return raw_error

def gaus(x, a, mu, sigma):
    return (a*1/ (sigma * np.sqrt(2*np.pi)) ) * np.exp( -0.5 * ((x-mu)/sigma)**2 )

#function used for minimizing. This means that, for a good fit, 
def chisqfunc(arguments, data):
    bins = data.get('args')[0][0]
    counts = data.get('args')[0][1]
    errors = data.get('args')[0][2]
    
    chisq = (counts - gaus(bins, *arguments))/errors
    return np.sum(chisq**2)

def plotting(bins, counts, counts_errors, result, popt, result_var):
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot()

    ax.errorbar(bins, counts, yerr = counts_errors, drawstyle='steps-mid', color='black')
    ax.plot(bins, gaus(bins, *result.x), label='Chi2 Minimization')
    ax.plot(bins, gaus(bins, *popt), label='LSQ Minimization')
    ax.grid()
   

    # #here we are building our text bubble to be plotted
    textstr = r'''$\chi2$ Minimization
    A = {:.0f} +/- {:.2f}
    $\mu$ = {:.2f} +/- {:.2f}
    $\sigma$ = {:.2f} +/- {:.2f}
    $\chi2$/dof = {:.1f}/{}'''.format(result.x[0], result_var[0],result.x[1], result_var[1],result.x[2], result_var[2], result.fun, len(counts)-3)



    at1 = AnchoredText(textstr, prop=dict(size=14), frameon=True, loc='upper right')
    at1.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax.add_artist(at1)
    
    (mu2, sig2) = norm.fit(data)
    
    y = norm.pdf(bins, mu2, sig2)
    y2 = minimize(chisqfunc, [A, mu2, sig2], args={'args':arguments})
    
    # ax.plot(bins, gaus(bins, *y2.x), 'b--.', label='Stats.norm values')
    
    
    
    ax.legend(loc='upper left', fontsize=15)
    plt.show()



#---------Globals----------#
freq_val = 5
pm = 'PM0'
pm_name = '..\PM0\Frequency'
pmt_filepath = '\Freq_{}kHz'.format(freq_val)
A = 450
mu = 309   
sigma=1.5


if __name__ == '__main__':
    
    #------Importing the data and changing units from C to nC
    df = pd.read_csv(pm_name+pmt_filepath+'\Data\charge_max_5_bi.csv')
    data = np.array(df.iloc[:, 0]) * 10**9
    
    #---------Performing the Histogram and changing the bin widths to the bins-----------#
    counts, bin_width = np.histogram(data, bins=25)
    bins = 0.5*(bin_width[1:] + bin_width[:-1])

    #-----------Setting inital guess, performing a lsq fit first, finding the appropriate errors, and then performing the chi2 Minimization--------------#
    ig = [A, mu, sigma]
    popt, pcov = curve_fit(gaus, bins, counts, p0=ig)

    counts_errors = error_bins(counts, gaus(bins, *popt), np.sqrt(counts))

    arguments = ([bins, counts, counts_errors], )
    result = minimize(fun=chisqfunc, x0=popt, args={'args':arguments})

    #-------------Finding the errors on the parameters--------------#
    # result_hess = result.hess_inv.matmat(np.eye(3))
    result_hess = result.hess_inv
    result_var = np.sqrt(np.diag(result_hess))


    print('LSQ chi2 value: {}'.format(np.sum( (counts - gaus(bins, *popt))**2/counts_errors**2)))

#-----------------Plotting the histogram and models----------------#
    plotting(bins, counts, counts_errors, result, popt, result_var)
