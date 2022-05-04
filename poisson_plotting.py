from math import factorial
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import scipy as sp
import math
from scipy.optimize import curve_fit
import pandas as pd

#df = pd.read_csv('..\{}\Data\{}.csv'.format("PMT0_waveform_1", " total_waveform_integrate_charge"))
df = pd.read_csv('..\{}\Data\{}.csv'.format("PMT0_waveform_5", "charge_max_5_bi"))
data = np.array(df[df.columns[0]].values.tolist())*10**9 #here the data is transformed into nC

fig, ax = plt.subplots(2, 1)


actual_count, actual_bins, patches = ax[0].hist(data, bins=70, histtype='step', range=(0, 20))
ax[0].set_title("Actual Data")
ax[0].set_xlabel('Charge [nC]')
ax[0].grid()
actual_bins = 0.5 * (actual_bins[1:] + actual_bins[:-1])



mu=4.5 #number of PE we expect. for toy model we can play around with this to our liking.

#num_of_events = np.arange(0, 20, 0.002)
num_of_events = np.random.poisson(lam=mu, size=10000)
#num_of_events = num_of_events//2
pe_charge = 1.602*10**-19

def poisson(k):
    return (10000*(mu)**k * math.exp(-mu))/math.factorial(k)

pe = np.arange(0, 10, 1) # here we are generating our number of PE, meaning we create an array that describes what number of PE we are looking at.
pe_list=[]
count=[]
for i in range(0, len(pe)):
    count.append(poisson(pe[i])) #here we are finding the number of counts (or expected events) for a certain number of PE. An example is passing through 1 PE will result in some count.
pe = pe*pe_charge*10**9 # puts the toy model into nano-coulombs



for i in range(1, len(pe)):
    random_gain = np.random.normal(2.031E10, 8.99E9, int(count[i])) # this will create a certain number of random numbers based on a normal distribution. The first number is the mean, the second if the std, and the final is the associated number of counts. An example is say we wanted a different gain associated with our one PE, then we pass in the counts for the PE and get out a different gain for that count. Each value within the count will have a different gain, and will be dictated by the random.normal() function.
    for j in range(0, len(random_gain)):
        pe_list.append(pe[i] * random_gain[j]) #here we are just multiplying the gain we found by the charge of PE. this will create a list of all measured charges within out toy model. 


def f(x, mu, Q1, sigma, N):
    #( (mu**1 * np.exp(-mu))/sp.special.factorial(1))
    #10**-8
    return N * ((( (mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma)) * np.exp((-1/(2*1)) * ((x-1*Q1)/sigma)**2)) + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma)) * np.exp((-1/(2*2)) * ((x-2*Q1)/sigma)**2)) + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma)) * np.exp((-1/(2*3)) * ((x-3*Q1)/sigma)**2)) + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma)) * np.exp((-1/(2*4)) * ((x-4*Q1)/sigma)**2)) + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma)) * np.exp((-1/(2*5)) * ((x-5*Q1)/sigma)**2)))

# plt.plot(num_of_events, y)
#counts, bins, patches = plt.hist(num_of_events, bins=10, range=(0, 10))
#plt.bar(pe, count)
counts, bins, patches = ax[1].hist(pe_list, bins=70, histtype='step', range=(0, 20))
bins = 0.5 * (bins[1:] + bins[:-1])
#, p0=[0.5, 50000000, 7500000, 10000]



props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

#, bounds=([0.5,65000000*pe_charge,9500000*pe_charge, 10**-10],[1.5,75000000*pe_charge,11500000*pe_charge, 1])
popt, pcov = curve_fit(f, bins, counts, p0=[0.1, 2.031E10*pe_charge*10**9, 8.99E9*pe_charge*10**9, 6243])
curve_fit_data = f(bins, *popt)
ax[1].plot(bins, curve_fit_data)
#ax[1].plot(bins, curve_fit_data)
print("Here are the returned values for the toy model")
print("mu={:.3e}, Q1={:.3e}, sigma={:.3e}, N={:.1f}\n".format(popt[0], popt[1]/(pe_charge*10**9), popt[2]/(pe_charge*10**9), popt[3]))
# plt.plot(x, y2)
ax[1].set_title("Toy Model")
ax[1].set_xlabel("Charge [nC]")
#n, bins, patches = plt.hist(num_of_events)
ax[1].grid()

#here we are curve fitting the actual data histogram
act_popt, act_pcov = curve_fit(f, actual_bins, actual_count, p0=[0.9, 2.031E10*pe_charge*10**9, 8.99E9*pe_charge*10**9, 6243])
ax[0].plot(actual_bins, f(actual_bins, *act_popt))

print("Here are the returned values for the actual data")
print("mu={:.3e}\nQ1={:.3e}\nsigma={:.3e}\nN={:.1f}".format(act_popt[0], act_popt[1]/(pe_charge*10**9), act_popt[2]/(pe_charge*10**9), act_popt[3]))


ax_text = "mu={:.3e}\nQ1={:.3e}\nsigma={:.3e}\nN={:.1f}".format(popt[0], popt[1]/(pe_charge*10**9), popt[2]/(pe_charge*10**9), popt[3])
ax[1].text(15, 300, ax_text, fontsize=14, verticalalignment='top', bbox=props)

ax0_text = "mu={:.3e}\nQ1={:.3e}\nsigma={:.3e}\nN={:.1f}".format(act_popt[0], act_popt[1]/(pe_charge*10**9), act_popt[2]/(pe_charge*10**9), act_popt[3])
ax[0].text(15, 500, ax0_text, fontsize=14, verticalalignment='top', bbox=props)




plt.show()
















