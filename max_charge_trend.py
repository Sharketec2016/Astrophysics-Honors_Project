'''
The purpose of this file is to take the extracted maximum charges from each pre-amplifier and plot a trend plot.


'''


import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import mean
from numpy.ma.core import std
import pandas as pd
import scipy as sp
from scipy.optimize import curve_fit
from scipy import stats
import math
from scipy import special
from numpy import asarray as ar, exp
from scipy.sparse import data


def linear(x, m, b):
    return m*x+b;



def open_file(pm_name, pm):
    pm_filepath = "..\{}\Frequency Change\{}".format(pm, pm_name)
    df = pd.read_csv('{}\Data\charge_max_5_bi.csv'.format(pm_filepath))
    data = np.array(df[df.columns[0]].values.tolist())
    return data



def plot(mean_charge, std_charge):
    
    frequency = [60, 70, 80, 90, 100]
   # plt.plot(frequency, mean_charge, 'o')
    #plt.show()
    popt, pcov = curve_fitting(mean_charge, frequency, std_charge)
    
    new_x = np.linspace(60, 100, 5)
    new_y = linear(new_x, *popt)
    
    chi2_x = np.linspace(60, 100, 5)
    chi2_y =linear(chi2_x, *popt)
    chi2_error = std_charge
    
   
    
    chi2_mean_charge = mean_charge
    chi2_star = np.sum(((chi2_mean_charge - chi2_y) / chi2_error) ** 2) / (len(chi2_mean_charge) - 2)
    
    #chi2, p = stats.chisquare(mean_charge, chi2_y, ddof=len(mean_charge)-2)
    
    plt.errorbar(frequency, mean_charge, yerr=std_charge, fmt='*')
    plt.plot(new_x, new_y, '--')
    textstr = 'b = {:.1f}+/-{:.2f}\nm = {:.1f}+/-{:.2f}\nchi2 = {:.4f}'.format(popt[1], pcov[1][1], popt[0], pcov[0][0], chi2_star)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    plt.text(85, 141, textstr, fontsize=14,
               verticalalignment='top', bbox=props)
    
    plt.xlabel("Frequency [kHz]")
    plt.ylabel("Charge [nC]")
    
    plt.title("Mean Max Charge for PM0")
    plt.grid()
    
    plt.show() 
    
    
    
def curve_fitting(mean_charge, frequency, std_charge):
    popt, pcov = curve_fit(linear, frequency, mean_charge, sigma=std_charge)
    return popt, pcov
        
    
    
def gain_calculation(mean_charge):
    #"\PM0\Frequency Change\PM0_60kHz_250mV\Data\ volts_max_ch_5_bi.csv"
    
    mean_input_voltage = []
    
    for i in range(6, 11):
        pm_filepath = "..\PM{}\Frequency Change\PM{}_{}0kHz_250mV".format("0", "0", i)
        df = pd.read_csv('{}\Data\{}'.format(pm_filepath, "int_max_volt.csv"))
        V_in_data = np.array(df[df.columns[0]].values.tolist())*10**-3 #here the units are Volts, from mV
        mean_input_voltage.append(V_in_data.mean())
    mean_input_voltage = np.array(mean_input_voltage)
    cap = 45 * 10**-12 #changed units to Farads [F]
    Q_in = cap * mean_input_voltage #[F]*[V] = [C]
    Q_in = Q_in * 10**9 
    
    
    
    Q_out = mean_charge
    gain = Q_out/Q_in #pC/pC
    
    gain_mean = gain.mean()
    gain_std = gain.std()
    
    return gain_mean, gain_std, gain
    




def main():
    
    mean_charge = []
    std_charge = []
    
    for i in range(6, 11):
        charge_data = open_file("\PM0_{}0kHz_250mV".format(i), 'PM0')*10**9
        mean_charge.append(charge_data.mean())
        std_charge.append(charge_data.std())
        
    plot(mean_charge=mean_charge, std_charge=std_charge)
    
    mean_gain, std_gain, gain = gain_calculation(mean_charge)
    print(mean_gain)
    
        
    
    
    
    
    
    
    
    
    
    
    
    
main()
