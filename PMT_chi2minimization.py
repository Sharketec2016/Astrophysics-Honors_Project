'''
The purpose of this file was to compare the minimization method from scipy.optimize. Scipy.optimize.minimize tries to minimize the least squares method, using a range of different
minimization methods. For the purpose of this work it was need to know whether or not a self implemented chi2 minimization would work. This file implements the same
minimization technique, that being BFGS, to minimize the objective function of LSQ and Chi2. 
'''



import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd
from matplotlib.offsetbox import AnchoredText
import fitting_functions_library as ff
from scipy.optimize import minimize

#----------Inital Guess Parameters---------#
gain = 0.8
std = 0.25
gain0 = 1
std0 = .1
alpha = 0.5
W = 0.5
N = 2500
N_ped = 8000


mu_df = pd.read_csv('Percentages of Pe peaks to total waveforms.csv')
mu_df = mu_df[mu_df.columns[5]].values.tolist()
mu_position = 2
temp_mu = mu_df[mu_position]

pois = np.random.poisson(temp_mu, 10000)
temp_counts, temp_bins = np.histogram(pois)

upper_limit=temp_bins[-1] + 3
lower_limit=0

#-------------------------------------------#
def lsq_method(func, bins, counts, p02, bnds):
    popt, pcov = curve_fit(func, bins, counts, p0=p02, bounds=bnds)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr


def grab_data(pmt_name):
    mu = mu_df[mu_position]
    df = pd.read_csv("./Charge Data/{}_{}.csv".format("total_integrate_charge", pmt_name))
    charge = np.array(df[df.columns[0]].values.tolist()) * 10**9
    return charge, mu


def par_update(parent_par, temp_par, name):
    temp_keys = list(temp_par.keys())
    
    if(name == 'Pedistal'):
        parent_par['name'] = name
        for k in parent_par:
            for i in range(len(temp_keys)):
                if(k == temp_keys[i]):
                    parent_par.get(k)[0] = temp_par.get(temp_keys[i])[0]
                    parent_par.get(k)[1] = temp_par.get(temp_keys[i])[1]
    elif(name == 'Noise'):
        parent_par['name'] = name
        for k in parent_par:
            for i in range(len(temp_keys)):
                if(k == temp_keys[i]):
                    parent_par.get(k)[0] = temp_par.get(temp_keys[i])[0]
                    parent_par.get(k)[1] = temp_par.get(temp_keys[i])[1]                
                    
                    

    return parent_par
    


def main():
    
    value = 0.1
    
    pmt_name = 'PMT0_2'
    charge, mu = grab_data(pmt_name)
    bin_width = np.arange(0, upper_limit+value, value)
    master_counts, master_bins= np.histogram(charge, bins=bin_width)
    master_bins = 0.5 * (master_bins[1:] + master_bins[:-1])
    
    plotting_counts, plotting_bins = master_counts, master_bins
    

    #--------------Here we are finding the cutting points for fitting the Pedistal, and the Noise--------------------#
    ped_upper = ff.ped_cut(master_bins, master_counts)[0]
    noise_lower = ff.noise_cut(master_bins, master_counts, mu)
    peaks_lower = ped_upper
    peaks_upper = noise_lower
    
    
    
    upper_bound = [np.inf, np.inf, np.inf, .8, 9, np.inf, 1, 20, np.inf]
    lower_bound = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001]
    initial_guess = [mu, std0, gain0, std, gain, N, W, alpha, N_ped]
    
    # popt, perr = lsq_method(ff.f, master_bins, master_counts, p02=initial_guess, bnds=(lower_bound, upper_bound))
    
    master_errors=ff.error_bins(master_counts, ff.f(master_bins, *initial_guess), np.sqrt(master_counts))
    
    pl_bins, pl_counts = master_bins, master_counts

    
    #----------------Use the data and plot the Pedistal------------------#
    master_counts, master_errors, ped_result, Q0_temp, ped_bins, ped_counts = ff.ped_handling(master_bins, master_counts, master_errors, initial_guess, ped_upper)
    ped_subtract_counts = master_counts
    
    #------------------Fitting and removing S_noise----------------------#
    master_counts, master_errors, noise_result, noise_bins, noise_counts = ff.noise_handling(master_bins, master_counts, master_errors, initial_guess, noise_lower, Q0_temp)
    noise_subtract_counts = master_bins
    #-----------------------S1 through S5 peaks being plotted--------------#
    master_counts, master_errors, peaks_result, peaks_bins, peaks_counts = ff.peaks_handling(master_bins, master_counts, master_errors, initial_guess, peaks_lower, peaks_upper)
    peaks_subtract_counts = master_counts
    #---------------------------------------------------------------------------#
    
    
    
    
    #--------------Computing the variance on the parameters--------------#
    
    #__________Pedistal__________#
    ped_hess = ped_result.hess_inv.matmat(np.eye(3))
    ped_var = np.sqrt(np.diag(ped_hess))

    #__________S-noise__________#
    noise_hess = noise_result.hess_inv.matmat(np.eye(3))
    noise_var = np.sqrt(np.diag(noise_hess))
    
    #__________S1-S5___________#
    peaks_hess = peaks_result.hess_inv.matmat(np.eye(5))
    peaks_var = np.sqrt(np.diag(peaks_hess))
    
    
    
    
    
    print('#---------Pedistal: sigma0, N, Q0, mu-------------#')
    print(ped_result)
    print('\nPedistal Error From Hessian: {}'.format(ped_var))
    print('\n#---------Noise: Q0, alpha, W, N, mu-------------#')
    print(noise_result)
    print('\nNoise Error From Hessian: {}'.format(noise_var))
    print('#---------Peaks: mu, Q0, sigma1, Q1, N-------------#')
    print(peaks_result)
    print('Peaks Error From Hessian: {}'.format(peaks_var))
    
    #---------------Populating the Parent Model for Plotting--------------------#
    arguments2 = [peaks_result.x[0], ped_result.x[0], ped_result.x[2], peaks_result.x[2], peaks_result.x[3], peaks_result.x[4], noise_result.x[2], noise_result.x[1], ped_result.x[1]]
    model = ff.f0(master_bins, *ped_result.x)+ff.s1_s5(master_bins, *peaks_result.x)+ff.S_noise(master_bins, *noise_result.x)
    data1 = ([master_bins, pl_counts, master_errors, model], )
    data = {'args': data1}
    chi2_value = ff.chisqfunc(arguments2, data)
    
    
    
    #----------------Grabbing the text bubble for Plotting----------------#
    #mu, sig0, q0, sig1, q1, N, W, alpha, chi2, dof, name
    textstr = ff.grab_text_bubble( [peaks_result.x[0], peaks_var[0]], [ped_result.x[0], ped_var[0]], [ped_result.x[2], ped_var[2]], [peaks_result.x[2], peaks_var[2]], [peaks_result.x[3], peaks_var[3]], [peaks_result.x[4], peaks_var[4]], [0, 0], [noise_result.x[1], noise_var[1]], [chi2_value], [len(master_bins) - 9], f"{pmt_name}")
    
    
    
    #---------------Plotting Histogram and Functions----------------#
    #________Figure and its attributes_______#
    fig = plt.figure(figsize=(20, 15))
    ax1 = fig.add_subplot()
    ax1.grid()
    #ax1.set_title(f'Data Set for {pmt_name}')
    ax1.set_title("Upper Limit for Plotting, and Noise cut off determined by Poisson Distribution")
    ax1.set_xlabel('Charge [nC]')
    ax1.set_ylabel('Counts')
    
    at1 = AnchoredText(textstr, prop=dict(size=14), frameon=True, loc='upper right')
    at1.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    ax1.add_artist(at1)
    
    #_________Plotting Histogram__________#
    ax1.errorbar(plotting_bins, plotting_counts, yerr=master_errors, drawstyle='steps-mid', color='black')
    #ax1.errorbar(master_bins, master_counts, yerr= master_errors, drawstyle='steps-mid', color='black')
    #_______Plotting Pedistal_____________#
    #ax1.plot(ped_bins, ff.f0(ped_bins, *ped_result.x), label='Pedistal')
    # #_______Plotting Noise_____________#
    #ax1.plot(noise_bins, ff.S_noise(noise_bins, *noise_result.x), label='Noise')
    # #_______Plotting Peaks_____________#
    #ax1.plot(peaks_bins, ff.s1_s5(peaks_bins, *peaks_result.x), label='Peaks')
    
    #master
    #mu, sigma0, Q0, sigma1, Q1, N, W, alpha, N_ped
    
    
   
    
    ax1.plot(master_bins, model)
    
    # print('#---------Pedistal: sigma0, N, Q0, mu-------------#')
    # print(ped_result)
    # print('\nPedistal Error From Hessian: {}'.format(ped_var))
    # print('\n#---------Noise: Q0, alpha, W, N, mu-------------#')
    # print(noise_result)
    # print('\nNoise Error From Hessian: {}'.format(noise_var))
    # print('#---------Peaks: mu, Q0, sigma1, Q1, N-------------#')
    # print(noise_result)
    # print('Peaks Error From Hessian: {}'.format(noise_var))
    
    
    ax1.set_xticks(np.arange(0, upper_limit, 1))
    plt.show()
    

#The order of the data will go as, Pedistal, Noise, and Peaks
#Columns: Chi2 Minimized Value, Errors from Hessian, Message, Number of Iterations, Convergence, Optimal Parameters
    
    
    
    






if __name__ == '__main__':
    main()
    








