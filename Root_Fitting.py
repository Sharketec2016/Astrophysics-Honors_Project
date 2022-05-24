#--------------------------------------------------------------#
'''
This file imports the electronic chain data. It takes the test pulse and output pulse maximum charges extracted from the integration from a previous file Read_Test.py. 
A histogram is created for the gain, and then fitted using a gaussian curve.
Parameters can be extracted from the fit if saving the values to a file is of interest.
The fitted gain histograms is saved as a pdf
The method for fitting the curve was done utilizing Minuit from ROOT, a C/C++ library created by CERN. For the purpose of this analysis, the python wrapper
of ROOT, PYROOT, was implemented to use the chi2 minimization method


'''
#--------------------------------------------------------------#
import ROOT as R
import numpy as np


#Email_Data/PM0/Frequency/PM0_60_kHz_250mV



def grab_input_data(pm_name, pm):

    pm_filepath = "../Data/{}/Frequency/{}".format(pm, pm_name)
    df = np.loadtxt('{}/total_charge.csv'.format(pm_filepath), delimiter=' ', dtype=float)
    df2 = np.loadtxt('{}/int_max_volt_fixed.csv'.format(pm_filepath), delimiter=' ', dtype=float)
    input_data = df
    test_data = df2*0.001#changing from mV to V
    return np.array(input_data), np.array(test_data), pm_filepath


#-------------------------------------------------------------------------------#
c1 = R.TCanvas('c1', 'Histogram', 5, 5, 800, 600)
R.gStyle.SetOptFit(1)
c1.SetHighLightColor(2);
c1.Range(143.029,-27.5625,146.9818,248.0625);
c1.SetFillColor(0);
c1.SetBorderMode(0);
c1.SetBorderSize(2);
c1.SetFrameBorderMode(0);
c1.SetFrameBorderMode(0);

pm_val = 0
freq_val = 100

pm = f'PM{pm_val}'
pm_name = f'PM{pm_val}_{freq_val}kHz_250mV' 

input_data, test_data, file_path = grab_input_data(pm_name, pm)

cap = 20e-12
test_data = cap*test_data
gain = input_data/test_data



bins = 30

h1 = R.TH1F(pm_name, 'Histogram', bins, gain.min()-5, gain.max()+5)
for i in range(len(gain)):
    h1.Fill(gain[i])


h2 = R.TH1F('input', 'Input_data', bins, input_data.min(), input_data.max())
for i in range(len(input_data)):
    h2.Fill(input_data[i])


h3 = R.TH1F('test', 'Test_data', bins, test_data.min(), test_data.max())
for i in range(len(test_data)):
    h3.Fill(test_data[i])



f1 = R.TF1('f1', '[0] * exp(-0.5*((x-[1])/[2]) * ((x-[1])/[2]))', gain.min()-5, gain.max()+5)

f1.SetParNames('Amplitude', 'Mean', 'Sigma', 'Flat')
f1.SetParameters(100, h1.GetMean(), h1.GetRMS())

h1.Fit('f1')
h1.SetXTitle('Gain')
#input()

'''
h2.Draw()
input()

h3.Draw()
input()
'''


pt = R.TPaveText(0.3978947,0.9345833,0.6021053,0.995,"blNDC");
pt.SetName("title");
pt.SetBorderSize(0);
pt.SetFillColor(0);
pt.SetFillStyle(0);
pt.SetTextFont(42);
pt_LaTex = pt.AddText("Histogram");
pt.Draw();
c1.Modified();
c1.cd();
c1.SetSelected(c1);

str = "{}/{}_Fixed.pdf".format(file_path, pm_name)
c1.Print(str)
#input()
