from math import factorial
from tkinter import W
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import scipy as sp
import math
from scipy.optimize import curve_fit
import pandas as pd
import unicodedata as ud


def h(x, sigma0, Q0, sigma1, Q1, N):
    
    mu = 0.1945
    return N*(  (1/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)  +
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
    

def f(x, sigma0, Q0, sigma1, Q1, N):
    
    mu = 0.1945
    return N*(  (1/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)  +
              + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1)/sigma1)**2)) 
              + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1)/sigma1)**2)) 
              + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1)/sigma1)**2)) 
              + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1)/sigma1)**2)) 
              + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1)/sigma1)**2))
              + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1)/sigma1)**2)) 
              + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1)/sigma1)**2))  
              + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1)/sigma1)**2)) 
              + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1)/sigma1)**2)) 
              + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1)/sigma1)**2)) )
        

pe_charge = 1.602*10**-19
pmt_name = "PMT0_waveform_3"
df = pd.read_csv("..\{}\Data_MB\{}.csv".format(pmt_name, " total_waveform_integrate_charge"))
charge = np.array(df[df.columns[0]].values.tolist()) * 10**9

upper_limit = 7

plt.figure(figsize=(10, 5))
counts, bins, ignore = plt.hist(charge, bins=140,range=(0, upper_limit))
bins = 0.5 * (bins[1:] + bins[:-1])

# #here are the parameters from the Monte Carlo Simulation.
mean= 0.1945
# gain = .4*10**10
# std = .01*10**9

# gain0 = 0.7*10**10
# std0 = .1*10**8

gain = 1.7
std = 0.09

gain0 = 1
std0 = .05

N = 10000


function_choice = 'f'

#here we need to limit the values on the parameters as to increase the stability of the model
unit_change = pe_charge*10**9

# if(function_choice == 'g'):
#     initial_guess = [mean, std0*unit_change, gain0*unit_change, std*unit_change, gain*unit_change, 5, 0.5]
# else:    
#     initial_guess = [std0, gain0, std, gain]
    
initial_guess = [std0, gain0, std, gain, N]
#here is the inital guess for the function f:

upper_bound = [np.inf, 1.5, np.inf, np.inf, np.inf]
lower_bound = [0, 0, 0, 1.5, 0]
#x, mu, sigma0, N, Q0, sigma1, Q1, N1
popt, pcov = curve_fit(h, bins, counts, p0=initial_guess)#, bounds=(lower_bound, upper_bound))

curve_fitted_data = h(bins, *popt)
plt.plot(bins, curve_fitted_data)
print(popt)
chi2 = np.sum((counts - curve_fitted_data)**2 / curve_fitted_data)
#chi2, p = sp.stats.chisquare(counts, curve_fitted_data, ddof=len(charge)-7)

#test_str = r'$\chi2$'

# if(function_choice == 'g'):
#     textstr = "$\mu$ = {:.3e}\nsimga0 = {:.3e}\nN = {:.3e}\nQ0 = {:.3e}\nsigma1 = {:.3e}\nQ1 = {:.3e}\n$\ alpha$ = {:.3e}\nW = {:.3e}\n$\chi2$/dof = {:.0f}/{}".format(popt[0], popt[1]/(pe_charge*10**9), popt[2], popt[3]/(pe_charge*10**9), popt[4]/(pe_charge*10**9), popt[5]/(pe_charge*10**9), popt[6], popt[7], chi2, 7)    
# else:
#     textstr = "simga0 = {:.3f}\nQ0 = {:.3f}\nsigma1 = {:.3f}\nQ1 = {:.3f}\n$\chi2$/dof = {:.0f}/{}".format(popt[0],popt[1], popt[2], popt[3], chi2, 7)

textstr = "simga0 = {:.3f}\nQ0 = {:.3f}\nsigma1 = {:.3f}\nQ1 = {:.3f}\nN = {:.3f}\n$\chi2$/dof = {:.0f}/{}".format(popt[0],popt[1], popt[2], popt[3], popt[4], chi2, 7)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

#place a text box in upper left in axes coords
plt.text(3, np.max(counts), textstr, fontsize=14,verticalalignment='top', bbox=props)


plt.xlabel('Charge [nC]')
plt.xticks([r for r in range(upper_limit+1)])
plt.title(pmt_name)
#plt.yscale('log')
plt.grid()
plt.show()













