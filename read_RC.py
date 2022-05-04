import numpy as np
import bias_RC
import charge_RC
import math
import pandas as pd

df = pd.DataFrame()
# this file reads all the raw data analyze it
# creating new files with the information

# define all the values needed for the functions and import all the two needed modules
pmt_path = '..\PMT0_waveform_5'
pmt_name = "PMT0_waveform_5"
time_spacing = 0.4
#number for window should follow: 3128 points per waveform. 3128/5 us = 625.6 points/us -> 625.6*3.5us = 2189.6 points
window_integration = 626  # size (in points) of the integration window for the charge. Found by finding the number of points for each unit of time. Then finding how many points
#for the total length of the pulse.
resistance = 50
files = 10000  # total has to be set by hand, for now
# number of traces saved as .csv files

# initialize some lists and arrays to store processed waveforms and integrated charge windows
v_list = []
v_array = []
ch_list = []
ch_array = []

# no longer needed, now grabbed from our picoscope generated .csv files
# time=np.arange(0,(points*time_spacing),time_spacing)
# creating the time array with the points and time spacing specifications

# open all the files needed to save the information analyze

charges = open(pmt_path + '\Data\charge_5_bi.csv', 'w')
# file for charges of the waveforms
charges_max = open(pmt_path +'\Data\charge_max_5_bi.csv', 'w')
# file for max charge of each waveform
base_line = open(pmt_path +'\Data\ base_line_5_bi.csv', 'w')
# file for baseline of each waveform
voltage_base = open(pmt_path +'\Data\ volt_base_5_bi.csv', 'w')
# file for processed waveforms
volts_max = open(pmt_path +'\Data\ volts_max_ch_5_bi.csv', 'w')
# file for the max amplitude of each waveform
time_ch = open(pmt_path + '\Data\ time_max_ch_5_bi.csv', 'w')
# file for the time when the max charge of each waveform occur
time_v = open(pmt_path + '\Data\ time_max_v_5_bi.csv', 'w')
# file for the time when the max amplitude of each waveform occur
max_all_time = open(pmt_path + '\Data\ max_ch_v_and_time.csv', 'w')
# file joining the baseline, the max amplitude, the time the max amplitude occur
int_volts = open(pmt_path +  '\Data\integrated_voltage.csv', 'w')
total_waveform_charge = open(pmt_path+'\Data\ total_waveform_integrate_charge.csv', 'w')
int_max_volt = open(pmt_path +  '\Data\int_max_volt.csv', 'w')

resistor_ch_max = open(pmt_path+'\Data\Resistor_charge_max.csv', 'w')
resistor_ch = open(pmt_path+'\Data\Resistor_charge.csv', 'w')
# the max charge and the time when the max charge occur
charges_maxx = []
base_linee = []
volts_maxx = []
time_chh = []
time_vv = []
range1 = files + 1
for file in range(1, range1):  # picoscope saves starting *_01
    print(file)
    n = str(file)

    #    start=waveform*points

    rfilename = '{0}\{1}_{2}.csv'.format(pmt_path, pmt_name, n.rjust(5, '0'))
    ini = np.genfromtxt(rfilename, delimiter=",", skip_header=3)
    # read in .csv
    #print(ini[:, 1])


   # volts = (ini[:, 1] * -1) * 2
    volts = ini[:, 1] * -1 * 93/50
    pre_chain_voltage= ini[:, 2]
    
    
    #here we want to integrate the current flowing over the resistor right by the LED, this will give an idea of how many photons should have been produced per pulse.
    resistor_voltage = ini[:, 3]

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

    base = bias_RC.bias(voltage=volts, slope=points)
    # print('Write baseline to file')
    # find the baseline
    # input list of amplitudes and the number of points in a waveform.
    base_line.write(str(base) + '\n')
    base_linee.append(base)
    # write baseline to file
    # print('Write baseline to file')

    volts_base, voltage_max, time_max_v = bias_RC.voltage_base_and_max(volts=volts, bias=base, time=time)
    # get an array of voltages without the baseline, get the maximum amplitude of that array and the time when that happens
    # input list of amplitudes, baseline calculated and the time spacing
    #volts_base = np.zeros(len(volts_base))
    v_list.append(volts_base)
    # print('called bias_RC')
    volts_max.write(str(voltage_max) + '\n')
    time_v.write(str(time_max_v) + '\n')
    volts_maxx.append(voltage_max)
    time_vv.append(time_max_v)
    # write everything to its respective file
    # print('Wrote everything to files')

    integrated_volts=[]
    integrated_max_volts=[]

    
    
    ch, charge_max, time_max_ch, Volts, V_max, total_wave_charge = charge_RC.charge(window=window_integration, time=time, volts=volts_base,
                                                   resistance=resistance, time_spacing=time_spacing, points=points, pre_chain_volts=pre_chain_voltage)
    # get an array of charge, the maximum charge in that array and the time when that happens
    # input the size of the window of integration, a list of the corrected amplitudes, the list of time, the time spacing, the value of resistance and the number of points in a waveform

    ch_list.append(ch)
    integrated_volts.append(Volts)


    #here we are passing in the voltage passing through the resistor to integrate it, and hopefully know the full charge passing through the LED. Window for points for integration was found to be 0.375 us, so window is 625.6*0.375 = 234.6 points
    ch_resistor, charge_max_resistor, time_max_ch_resistor, Volts_resistor, V_max_resistor, total_wave_charge_resistor = charge_RC.charge(window=235, time=time, volts=resistor_voltage, resistance=220, time_spacing=time_spacing, points=points, pre_chain_volts=pre_chain_voltage)


    resistor_ch_max.write(str(charge_max_resistor)+'\n')
    resistor_ch.write(str(ch_resistor)+'\n')
    total_waveform_charge.write(str(total_wave_charge) + '\n')

    charges_max.write(str(charge_max) + '\n')
    int_max_volt.write(str(pre_chain_voltage.max()) + '\n')
    time_ch.write(str(time_max_ch) + '\n')
    charges_maxx.append(charge_max)
    time_chh.append(time_max_ch)
    # write everything to its respective file
df['Number of NAN counts for each waveform file'] = count_of_nan
df.to_csv('NAN counts for each file.csv')

v_array = np.transpose(np.array(v_list))
np.savetxt(voltage_base, v_array, delimiter=",")

ch_array = np.transpose(np.array(ch_list))
np.savetxt(charges, ch_array, delimiter=",")

number_waveform = np.arange(0, files, 1)
np.savetxt(max_all_time, np.column_stack((number_waveform, base_linee, volts_maxx, time_vv, charges_maxx, time_chh)),
           delimiter=',',
           header='Waveform, Baseline (mV), Maximum Amplitude (mV), Time of max amplitude (ns), Max Charge (pC), Time of max charge (ns)')
# save the max amplitude, the time when the max amplitude happens, the max charge and the time when it happens and the wavefomr of each data all in one file

charges.close()
charges_max.close()
base_line.close()
voltage_base.close()
volts_max.close()
time_ch.close()
time_v.close()
max_all_time.close()
total_waveform_charge.close()

