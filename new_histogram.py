import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
from scipy import stats
import math
from scipy import special
from numpy import asarray as ar, exp

# The purpose of this code is to display and save the histograms for max charge and max amplitude


def gaus(x, a, mu, sigma):
    return a * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))


def pois(x, mu, a):
    return a * np.exp(-mu) * mu ** x / sp.special.factorial(x)

def grab_data(pm_name, pm):

    pm_filepath = "..\{}\Frequency Change\{}".format(pm, pm_name)
    df = pd.read_csv('{}\Data\charge_max_5_bi.csv'.format(pm_filepath))
    data = np.array(df[df.columns[0]].values.tolist())
    return data



def curve_fitting(x, y, fittype):
    if(fittype == 'Gaussian'):
        popt, pcov = curve_fit(gaus, x, y, p0=[y.max(), x.mean(), x.std()])
        return popt, pcov
    elif(fittype == 'Poisson'):
        popt, pcov = curve_fit(pois, x, y)
        return [popt, pcov]


def plot_making(x, raw_data, curve_fit, errors, fittype, stats, popt, charge_data):
    fig, ax = plt.subplots(1, 2, figsize=(12, 8))

#Here we are plotting the histogram and its curve fit. Along with the corresponding labels, and xticks.
    ax[0].bar(x, raw_data, width=0.03, edgecolor='black', yerr=errors)
    ax[0].plot(x, curve_fit, marker='o', linestyle='-', label='fit result', color='orange')
    x_axis_display = np.linspace(x.min(), x.max(), math.floor(len(x) / 2))
    ax[0].set_xticks(x_axis_display)
    ax[0].grid()
    ax[0].set_ylabel('Counts')
    ax[0].set_xlabel('Charge (nC)')
    ax[0].set_title(fittype)
    '''
    textstr = '\n'.join((
        r'$A=%.2f$' % (popt[0],),
        r'$\mathrm{mean}={:.2e}$' % (popt[1],),
        r'$\sigma={:.2e}f$' % (popt[2],),
        r'$\chi2=%.3f$' % (chi2)))
'''
    textstr = "A={:.2f}\nmean={:.2e}\nsigma={:.2e}\nchi2={:.2f}".format(popt[0], popt[1], popt[2], stats[0])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax[0].text(0.05, 0.95, textstr, transform=ax[0].transAxes, fontsize=14,
               verticalalignment='top', bbox=props)


#here we are plotting the residuals of the residuals of our raw_data and curve_fit
    residuals = raw_data - curve_fit

    ax[1].errorbar(x, residuals, yerr=errors, fmt='o')
    zer = np.zeros(len(x))
    ax[1].plot(x, zer, color='black')
    ax[1].set_title('Residuals plot')

#this here is whether or not we decide we would like tosave the image. The file path will need to be maually entered since its not passed through
    # plt.savefig(pmt_filepath+pm_name+'\Plots\ ' + type_of_curve+'_bins_'+str(k))
    # plt.savefig(pmt_filepath+pm_name+'\Plots\Adjusted Charge')

    #ax[0].plot(charge_data, gaus(charge_data, 150, charge_data.mean(), charge_data.std()), marker='o', linestyle='', color='orange')

    plt.show()


def error_bins(raw_data, fitted_data, raw_error):
    chi2_table = {0: [0.00, 1.29], 1: [0.27, 2.75], 2: [0.74, 4.25], 3: [1.10, 5.30], 4: [2.34, 6.78], 5: [2.75, 7.81],
                  6: [3.82, 9.28],
                  7: [4.25, 10.30], 8: [5.30, 11.32], 9: [6.33, 12.79], 10: [6.78, 13.81], 11: [7.81, 14.82],
                  12: [8.83, 16.29], 13: [9.28, 17.30],
                  14: [10.30, 18.32], 15: [11.32, 19.32], 16: [12.33, 20.80], 17: [12.79, 21.81], 18: [13.81, 22.82],
                  19: [14.82, 23.82], 20: [15.83, 25.30]}

    for i in range(len(raw_data)):
        if raw_data[i] < 10:
            error_band = chi2_table.get(raw_data[i])
            upper = abs(error_band[0] - raw_data[i])
            lower = abs(error_band[1] - raw_data[i])
            residual = fitted_data[i] - raw_data[i]
            if(residual >= 0):
                raw_error[i] = upper
            elif(residual < 0):
                raw_error[i] = lower
    return raw_error



def main():
    # here im going to grab the data from the csv. multiply by 10^9 to put us in units of pC
    charge_data = grab_data("\PM0_100kHz_250mV", 'PM0')*10**9

    '''
        # here im going backwards, and finding a list of the integrated. Im importing the list of the charges from each waveform, then im going to fnd the integrated voltages.
        df2 = pd.read_csv(pmt_filepath + pm_name + '\Data\int_max_volt.csv')
        # this is in mV, so we chage this to volts.
        voltage_data = np.array(df2[df2.columns[0]].values.tolist())
        # HERE im changing it from mV to V
        voltage_data = voltage_data * 0.001
        '''

    n, bins = np.histogram(charge_data, bins=15)

    bins = 0.5 * (bins[1:] + bins[:-1])
    domain = np.linspace(min(bins), max(bins), 15)
    # here im gathering the error on the counts  within the histogram. The error on the bin is the sqrt of the count, and anything less than 10
    # is set to an error of 1
    error = np.sqrt(n)




    popt, pcov = curve_fitting(domain, n, 'Gaussian')

    curve_fitted_data = gaus(domain, *popt)

    error_star = error_bins(n, curve_fitted_data, error)




    chi2_star = np.sum(((n - gaus(domain, *popt)) / error_star) ** 2) / (len(n) - 3)

    chi, p = sp.stats.chisquare(n, curve_fitted_data, ddof=len(n)-3)

    stats = [chi2_star, p]

    plot_making(domain, n, gaus(domain, *popt), error, 'Gaussian', stats, popt, charge_data)



    '''

        cap = 100 * 10**-12
        Q_in = cap * voltage_data

        Q_out = data

        gain = Q_out/Q_in

        gain_mean = gain.mean()
        gain_std = gain.std()
        gain_str = "{:e}".format(gain_mean)
        gain_std = "{:e}".format(gain_std)

        gain_document = open(pmt_filepath+pm_name+'\Data'+pm_name+'.txt', 'w')

        gain_document.write(gain_str + '\n')
        gain_document.write(gain_std + '\n')
        gain_document.close()



        gain_count, gain_bins = np.histogram(gain, bins=15)
        gain_domain = np.arange(0, 15)

        gain_bins = 0.5 * (gain_bins[1:] + gain_bins[:-1])


        gain_error = np.sqrt(gain_count)
        for i in range(len(gain_count)):
            if gain_count[i] < 10:
                gain_error[i] = 1

        gaus_gain_popt, gaus_gain_pcov = sp.optimize.curve_fit(gaus, gain_domain, gain_count)

        fig, ax = plt.subplots(1, 2, figsize=(12, 8))
        ax[0].bar(gain_bins, gain_count, width=10**2.5, edgecolor='black', yerr=gain_error)
        ax[0].plot(gain_bins, gaus(gain_domain, *gaus_gain_popt), marker='o', linestyle='-', color='orange')
        ax[0].grid()

        residuals = gain_count - gaus(gain_domain, *gaus_gain_popt)

        slope, intercept, r, p, se = sp.stats.linregress(gain_count, gaus(gain_domain, *gaus_gain_popt))
        rsq_line = slope*n2 + intercept
        #ax[1].plot(n2, rsq_line)

        ax[1].errorbar(gain_bins, residuals,  yerr=gain_error,  fmt='o')
        zer = np.zeros(len(gain_domain))
        ax[1].plot(gain_domain, zer, color='black')
        ax[1].set_title('Residuals plot')
        ax[1].grid()
        #plt.show()

        '''



main()

