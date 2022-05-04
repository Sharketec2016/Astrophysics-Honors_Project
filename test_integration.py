import numpy as np
from numpy.core.fromnumeric import mean 
import scipy.integrate as integrate
import cmath
import math
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.stats import chisquare


#file contains the code for the window integration to get the charge of each waveform

def charge (window, time, volts,resistance, time_spacing, points):
    #this function returns the list of charge related to the voltage list, also the maximum charge and the time when it happens
    #the inputs of this function are the window size for the integration (in ns), 
    #the time array for the integration, the amplitude array for the integration,
    #the resistance to calculate the charge, the time spacing and the number of points in the waveform

    window_points=int(window)
    
    #calculation of the window size with respect to the number of points in the waveform

    points_charge=int(points-window_points)
  
    #the amount of points that are going to be in the charge array

    ch=[]
    volt_int = []
    for s in range (points_charge):
        I=integrate.simps(volts[s:s+window_points],time[s:s+window_points])
       # print("I for iteration {}: {}".format(s, I))
        I2 = 0
        #integration of the waveform 
        V = I2
        #here im adding a 10^-6 because our time axis is in microseconds. I dont see anywhere where the code accounts for a time unit change
        q=(I/resistance) * 10**-6
        #calculation of the charge

        ch.append(q)
        volt_int.append(V)
        #append each individual charge to a list
    #print(ch)
    charges=np.array(ch)
    ch_max=charges.max()
    #get the max charge
    Volts = np.array(volt_int)
    V_max = Volts.max()
   
    time_max_ch=time[ch.index(ch_max)]
    #time when the max charge occur


    return ch, ch_max, time_max_ch, Volts, V_max



def neg_swap(volts):
    new_volts = []
    count=0
    for i in volts:
        if count > 550:
            if i<0:
                new_volts.append(0)
            else:
                new_volts.append(volts[count])
        else:
            new_volts.append(volts[count])
        count+=1

    return new_volts
    
    


def gain_calculation(pm_type):
    #"\PM0\Frequency Change\PM0_60kHz_250mV\Data\ volts_max_ch_5_bi.csv"
    df = pd.read_csv("..\{}\Data\charge_max_5_bi.csv".format(pm_type))
    output_charge = np.array(df[df.columns[0]].values.tolist())
    
    df2 = pd.read_csv("..\{}\Data\int_max_volt.csv".format(pm_type))
    return mean_gain, std_gain, mean_inputted_charge_list, mean_outputted_charge_list





def plotting(mean_gain, std_gain, mean_inputted_charge_list, mean_outputted_charge_list, data_type):
    
    
    if data_type == 'Frequency':
        x_axis = np.linspace(list(mean_gain.keys())[0], list(mean_gain.keys())[-1], len(list(mean_gain.keys())))
        popt, pcov = curve_fitting(list(mean_gain.keys()), list(mean_gain.values()), list(std_gain.values()))
        y_axis = linear(x_axis, *popt)
        
        chi2_star = np.sum(((list(mean_gain.values()) - y_axis) / list(std_gain.values())) ** 2) / (len(list(mean_gain.values())) - 2)
        
        
        
        plt.errorbar(list(mean_gain.keys()), list(mean_gain.values()),yerr=list(std_gain.values()), fmt='o')
        plt.plot(x_axis, y_axis)
        plt.xlabel("Frequency [kHz]")
        plt.ylabel("Average Gain")
        plt.title("Frequency Change")
        
        textstr = 'b = {:.1f}+/-{:.2f}\nm = {:.1f}+/-{:.2e}\nchi2 = {:.4f}'.format(popt[1], pcov[1][1], popt[0], pcov[0][0], chi2_star)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(85, 525000, textstr, fontsize=14,
                verticalalignment='top', bbox=props)
        
    elif data_type == 'Amplitude':
        x_axis = np.linspace(list(mean_gain.keys())[0], list(mean_gain.keys())[-1], len(list(mean_gain.keys())))
        popt, pcov = curve_fitting(list(mean_gain.keys()), list(mean_gain.values()), list(std_gain.values()))
        
        y_axis = linear(x_axis, *popt)
        #np.array(list(mean_gain.values())).mean()
        chi2_star = np.sum(((list(mean_gain.values()) - y_axis) / list(std_gain.values())) ** 2) / (len(list(mean_gain.values())) - 2)
        
        
        plt.errorbar(list(mean_gain.keys()), list(mean_gain.values()), yerr=list(std_gain.values()), fmt='o')
        plt.plot(x_axis, y_axis)
        #plt.hlines(y_axis, list(mean_gain.keys())[0], list(mean_gain.keys())[-1], colors='orange')
        plt.xlabel("Amplitude [mV]")
        plt.ylabel("Average Gain")
        plt.title("Amplitude Change")
        textstr = 'b = {:.1f}+/-{:.2f}\nm = {:.1f}+/-{:.2e}\nchi2 = {:.4f}\nMean Gain = {:.1f}'.format(popt[1], np.sqrt(pcov[1][1]), popt[0], np.sqrt(pcov[0][0]), chi2_star, np.array(list(mean_gain.values())).mean())
        #textstr = 'Mean Gain = {:.1f}\nb = {:.1f}+/-{:.2e}\nChi2 = {:.1f}'.format( np.array(list(mean_gain.values())).mean(), popt[0], np.sqrt(pcov[0][0]), chi2_star)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        plt.text(260, 9840, textstr, fontsize=14,
                verticalalignment='top', bbox=props)
        
        
        
    plt.grid()
    plt.show()
  


def linear(x, m, b):
    return m*x+b;

def B(x, b):
    return b

   
def curve_fitting(x, y, std):
    popt, pcov = curve_fit(linear, x, y, sigma=std)
   #popt, pcov = curve_fit(B, x, y, sigma=std)
    return popt, pcov
       



def main():
    
#here im going to be integrating the 1000 csv files myself, and compare marias integration against a trapazoid integration to see if it is consistent. 
   ''' 
    ini = np.genfromtxt('../PM0/Amplitude Change/PM0_100kHz_400mV/PM0_100kHz_400mV_0001.csv', delimiter=",", skip_header=3)
    volts = ini[:, 1] * 93/50
    pre_chain_voltage= ini[:, 2]
    # inverse polarity and correct for the impedance of the oscilloscope
    time = ini[:, 0]
    points = len(time)  # points in each waveform
    # print('grabbed the volts, time, and points')
    count_of_nan = 0
    nan_counts = []
    for i in range(0, len(volts)):
        if math.isnan(volts[i]):
            volts[i] = volts[i - 1]
            count_of_nan = count_of_nan + 1
        nan_counts.append(count_of_nan)

#here we are subtracting the background noise within the first 350 points. This is ok, since every pulse starts after the 350th element in the csv
    background_noise=volts[0:350].mean()
    volts = volts-background_noise



    #here we are obtaining the charge array, and the max_charge from this csv file.
    ch, ch_max, time_max_ch, Volts, V_max = charge(window=250, time=time, volts=volts, resistance=50, time_spacing=10, points=points)
    #here we are creating the new array from the original data that hase taken every value after the pulse, and set it to 0 if it were negative. Then we perform a trapezoid integration in the pulse. This integration should result in the same max charge as marias code. 
    new_volts_array = neg_swap(volts)
    std = np.array(new_volts_array).std()
    trapz = integrate.trapz(new_volts_array[0:1000],  dx=6.4)
    q=(trapz/50) * 10**-9
   
    
    
    fig, ax = plt.subplots(3, 1)
    ax[0].plot(ch)
    ax[0].set_ylabel("Charge")
    ax[0].set_title("Integrated Volts into Charge")
    
    ax[0].grid()
    
    ax[1].plot(new_volts_array)
    ax[1].set_ylabel("Volts")
    ax[1].set_title("New Volts Waveform After Changing Negative Into 0")
    ax[1].grid()
    #here i am passing  through the new voltage array into marias integration code. This should result in a value similar to the trapezoid integration
    ch, Q3, time_max_ch, Volts, V_max = charge(window=410, time=time, volts=new_volts_array, resistance=50, time_spacing=10, points=points)
    
    #here we are creating the table to be displayed
    val1 = ["Marias Code", "Trapaziod Integration", "New Array using Marias Integration"] #["{:X}".format(i) for i in range(2)] 
    val2 = ["{:02X}".format(10 * i) for i in range(1)] 
    val3 = [["{:.4e}".format(ch_max), "{:.4e}".format(q), "{:.4e}".format(Q3)] for r in range(1)] 
    ax[2].set_axis_off() 
    table = ax[2].table( 
        cellText = val3,  
        rowLabels = val2,  
        colLabels = val1, 
        rowColours =["palegreen"] * 10,  
        colColours =["palegreen"] * 10, 
        cellLoc ='center',  
        loc ='upper left',
        fontsize = 34)  
    table.auto_set_font_size(False)
    table.set_fontsize(12)    
    ax[2].set_title('Max Charge From Waveform', 
                fontweight ="bold") 
    
    plt.show()
    '''
    
    
   ''' 
    #here, after i have applied the new background subtraction on the waveforms, am plotting the mean gain against frequency. We should notice that the 
frequency=[60, 70, 80, 90, 100]
amplitude=[250, 300, 350, 400]
mean_gain_list = []
std_gain_list = []
print_delta_t=[]

Q_in_250=[]
Q_out_250=[]
print_delta_t_250=[]

#for i in range(6, 11):
for i in range(250, 450, 50):
    volt_df = pd.read_csv("../PM0/Amplitude Change/PM0_100kHz_{}mV/Data/int_max_volt.csv".format(i))
    charge_df = pd.read_csv("../PM0/Amplitude Change/PM0_100kHz_{}mV/Data/charge_max_5_bi.csv".format(i))

    volts_array = np.array(volt_df[volt_df.columns[0]].values.tolist())*10**-3
    charge_array = np.array(charge_df[charge_df.columns[0]].values.tolist())
    #new line
    delta_t = pd.read_csv("../PM0/Amplitude Change/PM0_100kHz_{}mV/Data/Delta_T_Risetime_TestPulse.csv".format(i))
    delta_t = np.array(delta_t[delta_t.columns[0]].values.tolist())
    Delta_T = delta_t
    delta_t = delta_t.mean()
    print_delta_t.append(delta_t)
    cap = 45 * 10**-12
    Q_in2 = cap*volts_array
    Q_in2 = Q_in2*10**9 #[nC]
    Q_out2 = charge_array*10**9
    
    if i == 350:
        Q_in_250 = Q_in2
        Q_out_250 = Q_out2
        print_delta_t_250 = Delta_T
    
    
    #new line
    gain = Q_out2/(Q_in2*delta_t)
    mean_gain2 = gain.mean()
    mean_gain_std = gain.std()
    mean_gain_list.append(mean_gain2)
    std_gain_list.append(mean_gain_std)

print(mean_gain_list)

plt.errorbar(amplitude, mean_gain_list, yerr=std_gain_list, fmt='*--')
plt.grid()
plt.xlabel("Amplitude")
plt.ylabel("Mean Gain")
plt.title("Mean Gain From Frequency")
plt.show()

count=0
for i in range(250, 450, 50):
    print("For a amplitude of {}, we have a delta T of: {}".format(i, print_delta_t[count]))
    count+=1

fig, axs = plt.subplots(3, 1)
axs[0].hist(Q_in2)
axs[0].grid()
axs[0].set_ylabel("Q_in2")
axs[1].hist(Q_out2)
axs[1].grid()
axs[1].set_ylabel("Q_out2")
axs[2].hist(Delta_T)
axs[2].grid()
axs[2].set_ylabel("Delta T")

plt.show()    

fig, axs = plt.subplots(3, 1)
axs[0].set_title("100kHz & 350 mV")
axs[0].hist(Q_in_250)
axs[0].grid()
axs[0].set_ylabel("Q_in_350")
axs[1].hist(Q_out_250)
axs[1].grid()
axs[1].set_ylabel("Q_out_350")


n, bins = np.histogram(print_delta_t_250)

bins = 0.5 * (bins[1:] + bins[:-1])
axs[2].bar(bins, n, width=0.01)
axs[2].grid()
axs[2].set_ylabel("Delta T 250")



plt.show()




'''



pm_type = 'PMT0_waveform_1'

mean_gain, std_gain, mean_inputted_charge_list, mean_outputted_charge_list = gain_calculation(pm_type=pm_type)

plotting(mean_gain, std_gain, mean_inputted_charge_list, mean_outputted_charge_list, data_type)


total_amplitude_mean_gain=[]
for key in mean_gain:
    print('{} has a gain of: {}'.format(key, mean_gain[key]))
    total_amplitude_mean_gain.append(mean_gain[key])
    
print("The mean gain for the entire PM setup for amplitude change is: {}".format(np.array(total_amplitude_mean_gain).mean()))
    
    
main()