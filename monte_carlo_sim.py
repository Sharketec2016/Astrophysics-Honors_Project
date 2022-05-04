import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare
from scipy.stats import binom
from scipy.stats import poisson
import scipy as sp






def f(x, mu, sigma0, N, Q0, sigma1, Q1):
    #( (mu**1 * np.exp(-mu))/sp.special.factorial(1))
    #10**-8
    return N*(  (1/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2 - mu)) 
              + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1)/sigma1)**2)) 
              + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1)/sigma1)**2)) 
              + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1)/sigma1)**2)) 
              + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1)/sigma1)**2)) 
              + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1)/sigma1)**2))  )






pe = [r for r in range(6)]
#this will create a list of PE looking like: [0, 1, 2, 3, 4, 5]
pe_charge=1.602*10**-19
#next we create the poisson distribution of the PE. We will say for this example that our expected value of PE is 1. Remember we cannot have a E(PE)=ration number, as 1.2 PE does not make any physical sense. 
mean=1
sigma = .3
waveforms=10000
pois_data = np.random.poisson(mean, waveforms)
data=np.ones(waveforms)
#now we will apply our gaussian to this poisson data. To do this, because out poisson stats will always be a integer value, we will pass each poisson value into our gaussian as its expected value. This is the process of convoluting our data.


for i in range(waveforms):
    if(pois_data[i] > 0):
        gaus_data = np.random.normal(pois_data[i], pois_data[i]*sigma, 1)
    else:
        gaus_data = np.random.normal(pois_data[i], sigma, 1)
    data[i] = gaus_data


#data = data * pe_charge * 10**9 #this line of code takes the histogram from before in PE, and changes it into charge on the x-axis. In particular it is in nC
plt.figure(figsize=(10, 5))
counts, bins, ignore = plt.hist(data, bins=90)
bins = 0.5 * (bins[1:] + bins[:-1])
#(mu, sigma0, N, Q0, sigma, Q1
popt, pcov = curve_fit(f, bins, counts, p0=[1, .1, 10000, 0, .3, 1], bounds=([0, 0, 0, 0, 0, 0], [5, .4, 100000, 0.001, .5, 5]))
#now that we have our optimized fitting parameters from our data to plug into our function f, we now create the curve fitted line


print("mu = {}".format(popt[0]))
print("sigma0 = {}".format(popt[1]))
print("N = {}".format(popt[2]))
print("Q0 = {}".format(popt[3]))
print("sigma1 = {}".format(popt[4]))
print("Q1 = {}".format(popt[5]))

new_domain = np.arange(bins[0], bins[-1], 0.0001)
fitted_data = f(new_domain, *popt)

#plt.xticks(pe)
plt.plot(new_domain, fitted_data)
textstr = "Mu = {:.3f}\nSimga0 = {:.3f}\nN = {:.3f}\nQ0 = {:.3f}\nSigma1 = {:.3f}\nQ1 = {:.3f}".format(popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
plt.text(7, 300, textstr, fontsize=14,verticalalignment='top', bbox=props)


plt.title("Set parameters: Mu={}, Sigma={}".format(mean, sigma))
plt.xlabel("Photo-electron")
plt.grid()
plt.xticks([r for r in range(8)])
#plt.savefig(".\Curve_fitting_plots\Mu={}_Sigma={}.jpeg".format(mean, sigma))
plt.show()
