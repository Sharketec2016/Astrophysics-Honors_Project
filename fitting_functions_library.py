import scipy as sp
import numpy as np
from scipy.optimize import minimize
from scipy.optimize import curve_fit



import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import poisson
import scipy as sp
import math
from scipy.optimize import curve_fit
import pandas as pd
from scipy.stats import chisquare
from matplotlib.offsetbox import AnchoredText

#--------------------Individual Peak, Overall Model, S_1 through S_8 Models-------------#
def f0(x, sigma0, N, Q0):
    return N*( (1)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2))

def f1(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2))

def f2(x, sigma1, N1, Q1, Q0, mu):
    return N1*( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 

def f3(x, sigma1, N1, Q1, Q0, mu):
    return N1*( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2))

def f4(x, sigma1, N1, Q1, Q0, mu):
    return N1*( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2))

def f5(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))

def f6(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2))

def f7(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1-Q0)/sigma1)**2))

def f8(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1-Q0)/sigma1)**2))

def f9(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1-Q0)/sigma1)**2))

def f10(x, sigma1, N1, Q1, Q0, mu):
    return N1*(((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1-Q0)/sigma1)**2))

def S_noise(x, Q0, alpha, N):
    return np.where(x<Q0, 0, N*alpha*np.exp(-alpha*(x)))

def f(x, mu, sigma0, Q0, sigma1, Q1, N, W, alpha, N_ped): 
    
    return np.where(x<=Q0, N*( ( ( (1)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2)) ) +   0*(alpha*W*np.exp(-alpha*(x-Q0) - mu))
            +   ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1)/sigma1)**2))
            + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1)/sigma1)**2)) 
            + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1)/sigma1)**2))  
            + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1)/sigma1)**2)) 
            + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1)/sigma1)**2)) 
            + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1)/sigma1)**2))), 
            
            
            N*( ( ( (1)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2)))   +   (alpha*W*np.exp(-alpha*(x-Q0) - mu))
            + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1)/sigma1)**2))
            + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1)/sigma1)**2)) 
            + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1)/sigma1)**2))  
            + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1)/sigma1)**2)) 
            + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1)/sigma1)**2)) 
            + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1)/sigma1)**2))) )
    
 
def s1_s5(x, mu, Q0, sigma1, Q1, N):
    return N*( ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2)) )
    

def error_bins(raw_data, fitted_data, raw_error):
    chi2_table = {0: [0.01, 1.29], 1: [0.27, 2.75], 2: [0.74, 4.25], 3: [1.10, 5.30], 4: [2.34, 6.78], 5: [2.75, 7.81],
                6: [3.82, 9.28],
                7: [4.25, 10.30], 8: [5.30, 11.32], 9: [6.33, 12.79], 10: [6.78, 13.81], 11: [7.81, 14.82],
                12: [8.83, 16.29], 13: [9.28, 17.30],
                14: [10.30, 18.32], 15: [11.32, 19.32], 16: [12.33, 20.80], 17: [12.79, 21.81], 18: [13.81, 22.82],
                19: [14.82, 23.82], 20: [15.83, 25.30]}

    for i in range(len(raw_data)):
        if raw_data[i] < 10:
            error_band = chi2_table.get(round(raw_data[i]))
            upper = abs(error_band[0] - raw_data[i])
            lower = abs(error_band[1] - raw_data[i])
            residual = fitted_data[i] - raw_data[i]
            if(residual >= 0):
                raw_error[i] = upper
            elif(residual < 0):
                raw_error[i] = lower
    return raw_error



def chisqfunc(arguments, data):
    bins = data.get('args')[0][0]
    counts = data.get('args')[0][1]
    errors = data.get('args')[0][2]
    model = data.get('args')[0][3]
    chisq = np.sum(((counts - model)**2/errors**2))
    return chisq



def chisqfunc_ped(arguments1, data):
    bins = data.get('args')[0][0]
    counts = data.get('args')[0][1]
    errors = data.get('args')[0][2]
    chisq1 = np.sum(((counts - f0(bins, *arguments1))**2/errors**2))
    return chisq1


def chisqfunc_noise(arguments2, data):
    bins = data.get('args')[0][0]
    counts = data.get('args')[0][1]
    errors = data.get('args')[0][2]
    chisq2 = np.sum(((counts - S_noise(bins, *arguments2))**2/errors**2))
    return chisq2

def chisqfunc_peaks(arguments3, data):
    bins = data.get('args')[0][0]
    counts = data.get('args')[0][1]
    errors = data.get('args')[0][2]
    chisq3 = np.sum(((counts - s1_s5(bins, *arguments3))**2/errors**2))
    return chisq3




def ped_handling(master_bins, master_counts, master_errors, ig, ped_upper):
    #mu, sigma0, Q0, sigma1, Q1, N, W, alpha, N_ped
    
    ped_upper_limit = ped_upper
    ped_lower_limit=0
    ped_bins = abs(master_bins[master_bins <= ped_upper_limit])
    ped_counts = abs(master_counts[master_bins <= ped_upper_limit])
    ped_errors = abs(master_errors[master_bins <= ped_upper_limit])
    
    initial_guess = [ig[1], ig[5], ig[2]]
    
    #these are the order of the needed input parameters: sigma0, N, Q0. Inital guesses from above were used, and the chi2 minimization was implemented
    ped_popt, ped_pcov = curve_fit(f0, ped_bins, ped_counts, p0=initial_guess, bounds=([0.001, 0.001, 0.001], [1, np.inf, 5]))
    ped_counts_errors = error_bins(ped_counts, f0(ped_bins, *ped_popt), np.sqrt(ped_counts))


    arguments = ([ped_bins, ped_counts, ped_counts_errors],)
    ped_result = minimize(chisqfunc_ped, ped_popt, bounds=((0.001, 1), (0.001, np.inf), (0.001, 1.5)), args={'args':arguments})
    #print(len(ped_counts))
   

    c = np.zeros( len(master_bins)-len(ped_bins))
    ped_bins_temp = np.concatenate((ped_bins, c))
    ped_counts_temp = np.concatenate((ped_counts, c))
    ped_counts_errors = error_bins(ped_counts_temp, f0(ped_bins_temp, *ped_result.x), np.sqrt(ped_counts_temp))
    Q0_temp = ped_result.x[2]
    
    master_counts = master_counts - f0(ped_bins_temp, *ped_result.x)
    master_errors = np.sqrt(master_errors**2 + ped_counts_errors**2)
    return master_counts, master_errors, ped_result, Q0_temp, ped_bins, ped_counts


def noise_handling(master_bins, master_counts, master_errors, ig, noise_lower, Q0):
    #mu, sigma0, Q0, sigma1, Q1, N, W, alpha, N_ped
    noise_lower_limit=noise_lower

    noise_bins = abs(master_bins[master_bins >= noise_lower_limit])
    noise_counts = abs(master_counts[master_bins >= noise_lower_limit])
    noise_errors = abs(master_errors[master_bins >= noise_lower_limit])


    initial_guess = [Q0, ig[7], ig[5]]
    #these are the order for the needed input parameters: Q0, alpha, N. Inital guesses from above were used. Chi2 minimization was implemented
    noise_popt, noise_pcov = curve_fit(S_noise, noise_bins, noise_counts, p0=initial_guess, bounds=([0.001, 0.001, 0.001], [5, 20, np.inf]))
    noise_counts_errors = error_bins(noise_counts, S_noise(noise_bins, *initial_guess), np.sqrt(noise_errors))
    arguments = ([noise_bins, noise_counts, noise_counts_errors],)
    noise_result = minimize(chisqfunc_noise, noise_popt, bounds=((0.001, 2), (0.001, 0.3), (0.001, 3000)), args={'args':arguments})


    c = np.zeros( len(master_bins)-len(noise_bins))
    noise_bins_temp = np.concatenate((noise_bins, c))
    noise_counts_temp = np.concatenate((noise_counts, c))
    noise_counts_errors = error_bins(noise_counts, S_noise(noise_bins, *noise_result.x), noise_counts_errors)
    master_counts = master_counts - S_noise(master_bins, *noise_result.x)

    master_errors[master_bins >= noise_lower_limit] = np.sqrt(master_errors[master_bins >= noise_lower_limit]**2 + noise_counts_errors**2)
    return master_counts, master_errors, noise_result, noise_bins, noise_counts


def peaks_handling(master_bins, master_counts, master_errors, ig, peaks_lower, peaks_upper):
    #mu,sigma0,Q0,sigma1, Q1, N, W, alpha
    peaks_lower_limit=peaks_lower
    peaks_upper_limit = peaks_upper

    peaks_bins = abs(master_bins[(master_bins > peaks_lower_limit) & (master_bins <= peaks_upper_limit)])
    peaks_counts = abs(master_counts[(master_bins > peaks_lower_limit) & (master_bins <= peaks_upper_limit)])
    peaks_errors = abs(master_errors[(master_bins > peaks_lower_limit) & (master_bins <= peaks_upper_limit)])
    #mu, Q0, sigma1, Q1, N.
    initial_guess = [ig[0], ig[2], ig[3], ig[4], ig[5]]
    #these are the order for the needed input parameters:  mu, Q0, sigma1, Q1, N. Inital guesses from above were used. Chi2 minimization was implemented
    print(initial_guess)
    peaks_popt, peaks_pcov = curve_fit(s1_s5, peaks_bins, peaks_counts, p0=initial_guess, bounds=([0.001, 0.001, 0.001, 0.001, 0.001], [6, 2, 1, 3, 2500]))
    
    peaks_counts_errors = error_bins(peaks_counts, s1_s5(peaks_bins, *peaks_popt), peaks_errors)
    arguments = ([peaks_bins, peaks_counts, peaks_counts_errors],)
    peaks_result = minimize(chisqfunc_peaks, peaks_popt, bounds=((0.001, 10), (0.001, 2), (0.001, 1), (0.001, 3), (0.001, 2500)), args={'args':arguments})
    
    
    return master_counts, master_errors, peaks_result, peaks_bins, peaks_counts


def ped_cut(bins, counts):
    cutting_pt = counts[(bins >= 0.8) & (bins <= 1.5)].min()
    cutting_pt = bins[np.where(counts == cutting_pt)[0]]
    return cutting_pt

def noise_cut(bins, counts, mu):

    # pois = np.random.poisson(mu, len(bins))

    # temp_counts, bins_temp= np.histogram(pois, density=True, bins=len(bins)//2)
    
    # return bins_temp[-1]
    
    cutting_pt = counts[bins > 2]
    cutting_pt = cutting_pt[cutting_pt < 25]
    cutting_pt = bins[np.where(counts == cutting_pt[0])[0]][0]
    return cutting_pt

def grab_text_bubble(mu, sig0, q0, sig1, q1, N, W, alpha, chi2, dof, name):
    #Note: all of the input parameters must be intered as either a tupple or a list that contains the variables value and its error to fill the text bubble. 
    textstr2 = r'''$\chi2$ Minimization of {}
    $\mu$ = {:.2f} +/- {:.2f}
    $\sigma_0$ = {:.2f} +/- {:.2f}
    $Q_0$ = {:.2f} +/- {:.2f}
    $\sigma_1$ = {:.2f} +/- {:.2f}
    $Q_1$ = {:.2f} +/- {:.2f}
    N = {:.0f} +/- {:.2f}
    W = {:.2f} +/- {:.2f}
    $\alpha$ = {:.2f} +/- {:.2f}
    $\chi2$/dof = ${:.1f}/{}$'''.format(name, mu[0], mu[1], sig0[0], sig0[1], q0[0], q0[1], sig1[0], sig1[1], q1[0], q1[1], N[0], N[1], W[0], W[1], alpha[0], alpha[1], chi2[0], dof[0])
    
    
    return textstr2



def main():
    print(f'Running __name__ = {__name__}')




if __name__ == '__main__':
    main()













