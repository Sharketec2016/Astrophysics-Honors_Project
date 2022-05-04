from audioop import mul
from math import factorial
from tkinter import W
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import scipy as sp
import math
from scipy.optimize import curve_fit
import pandas as pd
from scipy.stats import chisquare




def chi2_cal(charge, pois_bins):
    chi2_counts, chi2_bins = np.histogram(charge, bins=glob_bins,range=(lower_limit, upper_limit), density=True)
    chi2_bins = 0.5 * (chi2_bins[1:] + chi2_bins[:-1])
    
    try:
        chi2_popt, chi2_pcov = curve_fit(h, chi2_bins, chi2_counts)
    except RuntimeError:
        chi2 = np.nan
        p = np.nan
        return chi2, p
    curve_fitted_data = h(chi2_bins, *chi2_popt)
    normalized_model = curve_fitted_data * (np.sum(chi2_counts)/np.sum(curve_fitted_data))
    chi2, p = sp.stats.chisquare(chi2_counts, normalized_model, ddof=(len(chi2_counts)-7))

    
    man_chi2 = np.sum((chi2_counts - normalized_model)**2 / normalized_model)
    
    #print('chi2 :{:.3f} | manual chi2:{:8.3f}'.format(chi2, man_chi2))
    # bin_width = (chi2_bins.max() - chi2_bins.min()) / len(chi2_bins)
    # print(np.sum(chi2_counts * bin_width))
    # print(sp.integrate.simps(normalized_model, chi2_bins))


    return chi2/(len(chi2_counts)-7), p


def f(x, sigma0, Q0, sigma1, Q1, N, W, alpha): 
    return N*(  ( (1-W)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)  + (alpha*W*np.exp(-alpha*(x-Q0) - mu))
            + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))
            + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1-Q0)/sigma1)**2))  
            + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1-Q0)/sigma1)**2))) 
    
def h(x, sigma0, Q0, sigma1, Q1, N):
    return N*(  ( (1)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)
                + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))
                + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1-Q0)/sigma1)**2))  
                + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1-Q0)/sigma1)**2)) 
                + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1-Q0)/sigma1)**2)) )



mu_df = pd.read_csv('Percentages of Pe peaks to total waveforms.csv')
mu_df = mu_df[mu_df.columns[5]].values.tolist()
global mu 
mu = mu_df[6]
#4
pmt_name = "PMT1_waveform_4"
df = pd.read_csv("..\{}\Data_MB\{}.csv".format(pmt_name, " total_waveform_integrate_charge"))
charge = np.array(df[df.columns[0]].values.tolist()) * 10**9


#here i am setting the upper limit and lower limit for the fitting and drawing of the histogram. Also setting the global number of bins
upper_limit=7
lower_limit=0
glob_bins = 100

#defining the plot figure size
plt.figure(figsize=(15, 10))

#here we are plotting the histogram of the data
counts, bins, ignore = plt.hist(charge, bins=glob_bins,range=(lower_limit, upper_limit))
bins = 0.5 * (bins[1:] + bins[:-1])

#here are the inital guesses for the parameters
gain = 1.7
std = 0.9
gain0 = 1
std0 = .05
alpha = 0.2
W = 0.5

N = 10000
upper_bound = [np.inf, 50, np.inf, 9, np.inf]
lower_bound = [0, 0, 0, 0, 0]
initial_guess = [std0, gain0, std, gain, N]


#here we are splitting up the data
x_lessthan = []
x_greaterthan = []

y_lessthan=[]
y_greaterthan=[]

for i in range(len(bins)):
    if(bins[i] > 1):
        x_greaterthan.append(bins[i])
        y_greaterthan.append(counts[i])
    else:
        x_lessthan.append(bins[i])
        y_lessthan.append(counts[i])
    


#here we are curve fitting the different models to the different data sets. The arrays have been names appropriatly.
popt_lessthan, pcov_lessthan = curve_fit(f, x_lessthan, y_lessthan, p0=[std0, gain0, std, gain, N, alpha, W], bounds=([0, 0, 0, 0, 0, 0, 0], [np.inf, 1.5, np.inf, 9, np.inf, 1, 1]))
perr_lessthan = np.sqrt(np.diag(pcov_lessthan))

popt_greaterthan, pcov_greaterthan = curve_fit(h, x_greaterthan, y_greaterthan, p0=initial_guess, bounds=(lower_bound, upper_bound))
perr_greaterthan = np.sqrt(np.diag(pcov_greaterthan))



textstr = "{:30s}{}\nQ0 = {:.3f} +\- {:.3f} | Q0 = {:.3f} +\- {:.3f}\n$\sigma$_0 = {:.3f} +\- {:.3f} | $\sigma$_0 = {:.3f} +\- {:.3f}".format( "X < Q0", "| X >= Q0" ,popt_lessthan[1], perr_lessthan[1], popt_greaterthan[1], perr_greaterthan[1], popt_lessthan[0], perr_lessthan[0], popt_greaterthan[0], perr_greaterthan[0])


props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
plt.text(upper_limit-2, np.max(counts), textstr, fontsize=14,verticalalignment='top', bbox=props)


plt.title('Split Array | {}'.format(pmt_name))
plt.plot(x_greaterthan, h(x_greaterthan, *popt_greaterthan))
plt.plot(x_lessthan, f(x_lessthan, *popt_lessthan))

print(popt_greaterthan)
print(popt_lessthan)
#plt.yscale('log')
#plt.ylim(0, 10**2)
plt.grid()
#plt.savefig('.\Split_Array_plots\{}_bins_{}.jpeg'.format(pmt_name, glob_bins))
plt.show()



