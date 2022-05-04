import numpy as np 
import scipy.integrate as integrate
import cmath
import math

#file contains the code for the window integration to get the charge of each waveform

def charge (window, time, volts,resistance, time_spacing, points, pre_chain_volts):
    #this function returns the list of charge related to the voltage list, also the maximum charge and the time when it happens
    #the inputs of this function are the window size for the integration (in ns), 
    #the time array for the integration, the amplitude array for the integration,
    #the resistance to calculate the charge, the time spacing and the number of points in the waveform

    window_points=int(window)
    
    #calculation of the window size with respect to the number of points in the waveform

    points_charge=int(points-window_points)
  
    #the amount of points that are going to be in the charge array

    total_waveform_ch=[]

    I_tot = integrate.simps(volts[0:], time[0:])
    Q_tot = (I_tot/resistance) * 10**-6
    total_waveform_ch.append(Q_tot)
    total_waveform_ch = np.array(total_waveform_ch)
    total_waveform_ch = total_waveform_ch.max()
    ch=[]
    volt_int = []
    for s in range (points_charge):
        I=integrate.simps(volts[s:s+window_points],time[s:s+window_points])
        I2 = 0
        #integration of the waveform 
        V = I2
        #here im adding a 10^-6 because our time axis is in microseconds. I dont see anywhere where the code accounts for a time unit change
        q=(I/resistance) * 10**-6
        #calculation of the charge

        ch.append(q)
        volt_int.append(V)
        #append each individual charge to a list
    
    charges=np.array(ch)
    ch_max=charges.max()
    #get the max charge
    Volts = np.array(volt_int)
    V_max = Volts.max()
    
    
    
    
   
    time_max_ch=time[ch.index(ch_max)]
    #time when the max charge occur


    return ch, ch_max, time_max_ch, Volts, V_max, total_waveform_ch
