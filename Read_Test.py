'''
The purpose of this file was to input the csv waveform datasets for each test pulse setting, perform a background subtraction and then a total waveform integration
to extract the charge from each test pulse sent through the equipment. This extracted charge, from both the input and output test pulse, were saved into their
respective csv data files.
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import bias_RC
import charge_RC
import scipy.integrate as integrate

pmt_name = 'PM0_60kHz_250mV'
pmt_path = '..\PM0\Frequency Change\{}'.format(pmt_name)



resistance = 50
files = 1000 

total_waveform_charge = open(pmt_path+'\Data_Test\ total_charge.csv', 'w')
int_max_volt = open(pmt_path +  '\Data_Test\int_max_volt.csv', 'w')

for file in range(1, files+1):
    print(file)
    #here we are importing the data
    n = str(file)
    rfilename = '{0}\{1}_{2}.csv'.format(pmt_path, pmt_name, n.rjust(4, '0'))
    data = pd.read_csv(rfilename)
    test = np.array(data['Channel B'][3:]).astype(float)
    volts = np.array(data['Channel A'][3:]).astype(float) * 123/30
    time = np.array(data['Time'][3:]).astype(float)
    volts[np.isnan(volts)] = 0
    
    #here we are finding the baseline from the first 2us of data
    volts_base = bias_RC.bias(voltage=volts, slope=len(time))
    test_base = bias_RC.bias(voltage=test, slope=len(time))
    
    #here we are subtracting off the baseline from the data
    volts, voltage_max, time_max_v = bias_RC.voltage_base_and_max(volts=volts, bias=volts_base, time=time)

    test_voltage, pre_chain_max, time_max_pre_v = bias_RC.voltage_base_and_max(volts=test, bias=test_base, time=time)



    #here i am setting all points that are negative for the Blue curve to 0
    volts[469:][volts[469:] < 0] = 0


    #here we are integrating the entire waveform and finding the charge
    

    I_tot = integrate.simps(volts[0:], time[0:])
    Q_tot = (I_tot/resistance) * 10**-6 #here we are changing the units from us to seconds
    


    total_waveform_charge.write(str(Q_tot) + '\n')
    int_max_volt.write(str(test_voltage.max()) + '\n')




