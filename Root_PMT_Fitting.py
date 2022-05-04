# --------------------------------------------------------------#
'''
This file is used for importing, fitting, and presenting a fit of the PMT data. 
This file uses ROOT for fitting and plotting models and histograms of the data.

'''
# --------------------------------------------------------------#


from ROOT import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def grab_input_data(pm):
    pm_filepath = "../Charge PMT Data".format(pm)
    df = pd.read_csv('{}/total_integrate_charge_{}.csv'.format(pm_filepath, pm))

    return np.array(df.iloc[:, 0]) * 10 ** 9


# N, sigma0, Q0
def Sped(x, par):
    fn = par[0] / (TMath.Sqrt(2 * TMath.Pi()) * par[1]) * TMath.Exp(-0.5 * TMath.Sq((x[0] - par[2]) / par[1]))

    return fn


# "N", "Q0", "alpha" 'mu'
def Snoise(x, par):
    fn = par[0] * par[2] * TMath.Exp(-par[2] * (x[0] - par[1]) - par[3])
    return fn


def SN(x, par):  # N, Q0, Q1, mu, sigma1
    fn = par[0] * (par[3] * TMath.Exp(-par[3]) / (TMath.Sqrt(2 * TMath.Pi()) * par[4]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[2] - par[1]) / par[4]))  # S_1
                   + TMath.Power(par[3], 2) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(2) * 2 * par[4]) * TMath.Exp(
                -0.5 / 2 * TMath.Sq((x[0] - 2 * par[2] - par[1]) / par[4]))  # S_2
                   + TMath.Power(par[3], 3) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(3) * 3 * par[4]) * TMath.Exp(
                -0.5 / 3 * TMath.Sq((x[0] - 3 * par[2] - par[1]) / par[4]))  # S_3
                   + TMath.Power(par[3], 4) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(4) * 4 * par[4]) * TMath.Exp(
                -0.5 / 4 * TMath.Sq((x[0] - 4 * par[2] - par[1]) / par[4]))  # S_4
                   + TMath.Power(par[3], 5) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(5) * 5 * par[4]) * TMath.Exp(
                -0.5 / 5 * TMath.Sq((x[0] - 5 * par[2] - par[1]) / par[4]))  # S_5
                   + TMath.Power(par[3], 6) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(6) * 6 * par[4]) * TMath.Exp(
                -0.5 / 6 * TMath.Sq((x[0] - 6 * par[2] - par[1]) / par[4]))  # S_6
                   + TMath.Power(par[3], 7) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(7) * 7 * par[4]) * TMath.Exp(
                -0.5 / 7 * TMath.Sq((x[0] - 7 * par[2] - par[1]) / par[4]))  # S_7
                   + TMath.Power(par[3], 8) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(8) * 8 * par[4]) * TMath.Exp(
                -0.5 / 8 * TMath.Sq((x[0] - 8 * par[2] - par[1]) / par[4]))  # S_8
                   + TMath.Power(par[3], 9) * TMath.Exp(-par[3]) / (
                           TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(9) * 9 * par[4]) * TMath.Exp(
                -0.5 / 9 * TMath.Sq((x[0] - 9 * par[2] - par[1]) / par[4])))  # S_9
    return fn


# mu_peaks, sigma0, Q0, sigma1, Q1, Npeaks, alpha, Nped, Nnoise
def SParent(x, par):
    fn1 = (par[7]) / (TMath.Sqrt(2 * TMath.Pi()) * par[1]) * TMath.Exp(-0.5 * TMath.Sq((x[0] - par[2]) / par[1]))
    + par[5] * (par[0] * TMath.Exp(-par[0]) / (TMath.Sqrt(2 * TMath.Pi()) * par[3]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[4] - par[2]) / par[3]))  # S_1
                + TMath.Power(par[0], 2) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(2) * 2 * par[3]) * TMath.Exp(
                -0.5 / 2 * TMath.Sq((x[0] - 2 * par[4] - par[2]) / par[3]))  # S_2
                + TMath.Power(par[0], 3) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(3) * 3 * par[3]) * TMath.Exp(
                -0.5 / 3 * TMath.Sq((x[0] - 3 * par[4] - par[2]) / par[3]))  # S_3
                + TMath.Power(par[0], 4) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(4) * 4 * par[3]) * TMath.Exp(
                -0.5 / 4 * TMath.Sq((x[0] - 4 * par[4] - par[2]) / par[3]))  # S_4
                + TMath.Power(par[0], 5) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(5) * 5 * par[3]) * TMath.Exp(
                -0.5 / 5 * TMath.Sq((x[0] - 5 * par[4] - par[2]) / par[3]))  # S_5
                + TMath.Power(par[0], 6) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(6) * 6 * par[3]) * TMath.Exp(
                -0.5 / 6 * TMath.Sq((x[0] - 6 * par[4] - par[2]) / par[3]))  # S_6
                + TMath.Power(par[0], 7) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(7) * 7 * par[3]) * TMath.Exp(
                -0.5 / 7 * TMath.Sq((x[0] - 7 * par[4] - par[2]) / par[3]))  # S_7
                + TMath.Power(par[0], 8) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(8) * 8 * par[3]) * TMath.Exp(
                -0.5 / 8 * TMath.Sq((x[0] - 8 * par[4] - par[2]) / par[3]))  # S_8
                + TMath.Power(par[0], 9) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(9) * 9 * par[3]) * TMath.Exp(
                -0.5 / 9 * TMath.Sq((x[0] - 9 * par[4] - par[2]) / par[3])))  # S_9

    fn2 = par[7] / (TMath.Sqrt(2 * TMath.Pi()) * par[1]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[2]) / par[1]))  # Spedistal
    + par[8] * par[6] * TMath.Exp(-par[6] * (x[0] - par[2]) - par[0])  # Snoise
    + par[5] * (par[0] * TMath.Exp(-par[0]) / (TMath.Sqrt(2 * TMath.Pi()) * par[3]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[4] - par[2]) / par[3]))  # S_1
                + TMath.Power(par[0], 2) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(2) * 2 * par[3]) * TMath.Exp(
                -0.5 / 2 * TMath.Sq((x[0] - 2 * par[4] - par[2]) / par[3]))  # S_2
                + TMath.Power(par[0], 3) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(3) * 3 * par[3]) * TMath.Exp(
                -0.5 / 3 * TMath.Sq((x[0] - 3 * par[4] - par[2]) / par[3]))  # S_3
                + TMath.Power(par[0], 4) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(4) * 4 * par[3]) * TMath.Exp(
                -0.5 / 4 * TMath.Sq((x[0] - 4 * par[4] - par[2]) / par[3]))  # S_4
                + TMath.Power(par[0], 5) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(5) * 5 * par[3]) * TMath.Exp(
                -0.5 / 5 * TMath.Sq((x[0] - 5 * par[4] - par[2]) / par[3]))  # S_5
                + TMath.Power(par[0], 6) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(6) * 6 * par[3]) * TMath.Exp(
                -0.5 / 6 * TMath.Sq((x[0] - 6 * par[4] - par[2]) / par[3]))  # S_6
                + TMath.Power(par[0], 7) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(7) * 7 * par[3]) * TMath.Exp(
                -0.5 / 7 * TMath.Sq((x[0] - 7 * par[4] - par[2]) / par[3]))  # S_7
                + TMath.Power(par[0], 8) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(8) * 8 * par[3]) * TMath.Exp(
                -0.5 / 8 * TMath.Sq((x[0] - 8 * par[4] - par[2]) / par[3]))  # S_8
                + TMath.Power(par[0], 9) * TMath.Exp(-par[0]) / (
                            TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(9) * 9 * par[3]) * TMath.Exp(
                -0.5 / 9 * TMath.Sq((x[0] - 9 * par[4] - par[2]) / par[3])))  # S_9


def Fitting_Sped_SN(x, par):
    fn1 = (par[7]) / (TMath.Sqrt(2 * TMath.Pi()) * par[1]) * TMath.Exp(-0.5 * TMath.Sq((x[0] - par[2]) / par[1]))
    + par[5] * (par[0] * TMath.Exp(-par[0]) / (TMath.Sqrt(2 * TMath.Pi()) * par[3]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[4] - par[2]) / par[3]))  # S_1
                + TMath.Power(par[0], 2) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(2) * 2 * par[3]) * TMath.Exp(
                -0.5 / 2 * TMath.Sq((x[0] - 2 * par[4] - par[2]) / par[3]))  # S_2
                + TMath.Power(par[0], 3) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(3) * 3 * par[3]) * TMath.Exp(
                -0.5 / 3 * TMath.Sq((x[0] - 3 * par[4] - par[2]) / par[3]))  # S_3
                + TMath.Power(par[0], 4) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(4) * 4 * par[3]) * TMath.Exp(
                -0.5 / 4 * TMath.Sq((x[0] - 4 * par[4] - par[2]) / par[3]))  # S_4
                + TMath.Power(par[0], 5) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(5) * 5 * par[3]) * TMath.Exp(
                -0.5 / 5 * TMath.Sq((x[0] - 5 * par[4] - par[2]) / par[3]))  # S_5
                + TMath.Power(par[0], 6) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(6) * 6 * par[3]) * TMath.Exp(
                -0.5 / 6 * TMath.Sq((x[0] - 6 * par[4] - par[2]) / par[3]))  # S_6
                + TMath.Power(par[0], 7) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(7) * 7 * par[3]) * TMath.Exp(
                -0.5 / 7 * TMath.Sq((x[0] - 7 * par[4] - par[2]) / par[3]))  # S_7
                + TMath.Power(par[0], 8) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(8) * 8 * par[3]) * TMath.Exp(
                -0.5 / 8 * TMath.Sq((x[0] - 8 * par[4] - par[2]) / par[3]))  # S_8
                + TMath.Power(par[0], 9) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(9) * 9 * par[3]) * TMath.Exp(
                -0.5 / 9 * TMath.Sq((x[0] - 9 * par[4] - par[2]) / par[3])))  # S_9
    return fn1


def Fitting_Snoise_SN(x, par):
    fn2 = par[7] / (TMath.Sqrt(2 * TMath.Pi()) * par[1]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[2]) / par[1]))  # Spedistal
    + par[8] * par[6] * TMath.Exp(-par[6] * (x[0] - par[2]) - par[0])  # Snoise
    + par[5] * (par[0] * TMath.Exp(-par[0]) / (TMath.Sqrt(2 * TMath.Pi()) * par[3]) * TMath.Exp(
        -0.5 * TMath.Sq((x[0] - par[4] - par[2]) / par[3]))  # S_1
                + TMath.Power(par[0], 2) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(2) * 2 * par[3]) * TMath.Exp(
                -0.5 / 2 * TMath.Sq((x[0] - 2 * par[4] - par[2]) / par[3]))  # S_2
                + TMath.Power(par[0], 3) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(3) * 3 * par[3]) * TMath.Exp(
                -0.5 / 3 * TMath.Sq((x[0] - 3 * par[4] - par[2]) / par[3]))  # S_3
                + TMath.Power(par[0], 4) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(4) * 4 * par[3]) * TMath.Exp(
                -0.5 / 4 * TMath.Sq((x[0] - 4 * par[4] - par[2]) / par[3]))  # S_4
                + TMath.Power(par[0], 5) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(5) * 5 * par[3]) * TMath.Exp(
                -0.5 / 5 * TMath.Sq((x[0] - 5 * par[4] - par[2]) / par[3]))  # S_5
                + TMath.Power(par[0], 6) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(6) * 6 * par[3]) * TMath.Exp(
                -0.5 / 6 * TMath.Sq((x[0] - 6 * par[4] - par[2]) / par[3]))  # S_6
                + TMath.Power(par[0], 7) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(7) * 7 * par[3]) * TMath.Exp(
                -0.5 / 7 * TMath.Sq((x[0] - 7 * par[4] - par[2]) / par[3]))  # S_7
                + TMath.Power(par[0], 8) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(8) * 8 * par[3]) * TMath.Exp(
                -0.5 / 8 * TMath.Sq((x[0] - 8 * par[4] - par[2]) / par[3]))  # S_8
                + TMath.Power(par[0], 9) * TMath.Exp(-par[0]) / (
                        TMath.Sqrt(2 * TMath.Pi()) * TMath.Factorial(9) * 9 * par[3]) * TMath.Exp(
                -0.5 / 9 * TMath.Sq((x[0] - 9 * par[4] - par[2]) / par[3])))  # S_9
    return fn2


def Stotal(x,par):  # 'mu', 'sigma0', 'Q0', 'sigma1', 'Q1', 'N', 'alpha', 'N_ped', 'Nnoise' 'Q0noise', 'muNoise', 'Q0peaks

    # N, sigma0, Q0
    sped_par = [par[7], par[1], par[2]]

    ###"N", "Q0", "alpha" 'mu'
    snoise_par = [par[8], par[9], par[6], par[10]]

    # N, Q0, Q1, mu, sigma1
    sN_par = [par[5], par[2], par[4], par[0], par[3]]

    fn1 = Sped(x, sped_par) + Snoise(x, snoise_par) + SN(x, sN_par)
    fn2 = Sped(x, sped_par) + SN(x, sN_par)


    # if (x[0] > par[2]):
    #     return fn1
    # else:
    #     return fn2
    #return np.where(x[0] > par[2], fn1, fn2)
    return np.where(x[0] > par[2], np.where(x[0] < 6, Sped(x, sped_par) + SN(x, sN_par), Snoise(x, snoise_par) + SN(x, sN_par)), Sped(x, sped_par))


def Fit_Sped(ped_hist, n, ped_upperlim, ig):
    fit2 = TF1('fit2', Sped, -1, ped_upperlim, 3)
    fit2.SetParNames("N", "sigma0", "Q0")
    fit2.SetParameters(ig[0], ig[1], ig[2])
    fit2.SetParLimits(0, 1, n * 10)  # N
    fit2.SetParLimits(1, 0.001, 1)  # sigma0
    fit2.SetParLimits(2, 0.001, 1.5)  # Q0
    ped_hist.Fit("fit2", 'R')
    ped_hist.SetFillColor(7)

    return ped_hist, fit2


def Fit_Snoise(noise_hist, noise_lowerlim, ig, upper_fitlim):
    fit3 = TF1("fit3", Snoise, noise_lowerlim, upper_fitlim, 4)
    fit3.SetParNames("N", "Q0", "alpha", "mu")
    fit3.SetParameters(ig[0], ig[1], ig[2], ig[3])
    fit3.SetParLimits(0, 0, 10000)  # N
    fit3.SetParLimits(1, 0.001, 5)  # Q0
    # fit3.SetParLimits(1, ig[1], ig[1]);		#Q0

    fit3.SetParLimits(2, 0., 1)  # alpha
    fit3.SetParLimits(3, 0.001, 15)  # mu
    noise_hist.Fit("fit3", "R")
    print(f'Noise Fit Function At 8: {fit3.Eval(8)}')
    return noise_hist, fit3


# mu, Q0, sigma1, Q1, N
# bounds=((0.001, 10), (0.001, 2), (0.001, 1), (0.001, 3), (0.001, 2500))

def Fit_SN(peaks_hist, peaks_lowerlim, peaks_upperlim, ig):
    fit4 = TF1("fit4", SN, peaks_lowerlim, peaks_upperlim, 5)
    fit4.SetParNames("N", "Q0", "Q1", "mu", "sigma1")
    fit4.SetParameters(ig[0], ig[1], ig[2], ig[3], ig[4])
    fit4.SetParLimits(0, 0.5, 10000)  # N
    fit4.SetParLimits(1, 0.001, 5)  # Q0
    fit4.SetParLimits(3, 0.001, 10)  # mu
    fit4.SetParLimits(2, 0.001, 5)  # Q1
    fit4.SetParLimits(4, 0.001, 3)  # sigma1
    peaks_hist.Fit("fit4", "R")

    return peaks_hist, fit4


def Stotal_python(x, par):  # 'mu', 'sigma0', 'Q0', 'sigma1', 'Q1', 'N', 'alpha', 'N_ped', 'Nnoise'

    # N, sigma0, Q0
    sped_par = [par[7], par[1], par[2]]

    ###"N", "Q0", "alpha" 'mu'
    snoise_par = [par[8], par[2], par[6], par[0]]

    # N, Q0, Q1, mu, sigma1
    sN_par = [par[5], par[2], par[4], par[0], par[3]]
    #
    fn1 = Sped(x, sped_par) + Snoise(x, snoise_par) + SN(x, sN_par)
    fn2 = Sped(x, sped_par) + SN(x, sN_par)
    # if (x[0] > par[2]):
    #     return fn1
    # else:
    #     return fn2

    # fn1 = Fitting_Snoise_SN(x, par)
    # fn2 = Fitting_Sped_SN(x, par)

    return np.where(x > par[2], fn1, fn2)


# -------------------------------------------------------------------------------#

def main():
    # here we are creating the canvas for plotting
    c1 = TCanvas('c1', 'Histogram', 5, 5, 800, 600)
    c1.SetHighLightColor(2)
    c1.Range(143.029, -27.5625, 146.9818, 248.0625)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetBorderSize(2)
    c1.SetFrameBorderMode(0)
    c1.SetFrameBorderMode(0)

    # here i am grabbing the data
    pm = 'PMT1_1'
    wave = 3
    data = grab_input_data(pm)

    # here i am creating and populating the histogram
    bins = 90
    n = 0
    noise_lowerlim = 5

    ped_upperlim = 1.2

    upper_fitlim = 15
    mu_ig = 2.8

    parent_hist = TH1F('hist', pm, bins, -1, upper_fitlim)
    for i in range(len(data)):
        parent_hist.Fill(data[i])
        n += 1

    # here i am attempting to fit the pedistal and after fitting set the values to the variables
    # "N", "sigma0", "Q0"
    ped_hist = parent_hist.Clone()

    ped_ig = [n / 2, 0.5, 0]
    ped_hist, ped_fit = Fit_Sped(ped_hist, n, ped_upperlim, ped_ig)

    # here draw the pedistal peak
    ped_hist.Draw()
    input()

    N = ped_fit.GetParameter(0)
    sigma0 = ped_fit.GetParameter(1)
    Q0 = ped_fit.GetParameter(2)

    Nerr = ped_fit.GetParError(0)
    sigma0err = ped_fit.GetParError(1)
    Q0err = ped_fit.GetParError(2)

    # here i am fitting the noise term after i have subtracted off the pedistal
    noise_hist = parent_hist.Clone()
    ped_subtraction_fit = TF1('ped_sub', Sped, -1, ped_upperlim, 3)
    ped_subtraction_fit.SetParameters(N, sigma0, Q0)

    noise_hist.Add(ped_subtraction_fit, -1)

    # mu values from excel sheet will need to be manually entered
    # "N", "Q0", "alpha" "mu"
    noise_ig = [150, Q0, 0.5, 0.1]
    noise_hist, noise_fit = Fit_Snoise(noise_hist, noise_lowerlim, noise_ig, upper_fitlim)

    # here draw the noise fit
    noise_hist.Draw()
    input()

    N3 = noise_fit.GetParameter(0)
    alpha = noise_fit.GetParameter(3)
    Q0_noise = noise_fit.GetParameter(2)
    mu_noise = noise_fit.GetParameter(4)

    N3err = noise_fit.GetParError(0)
    alphaerr = noise_fit.GetParError(3)
    Q0_noise_err = noise_fit.GetParError(2)
    mu_noise_err = noise_fit.GetParError(4)

    # here i am creating a noise function to them subtract off of the histogram. After this i will fit the remaining SN peaks with a function
    peaks_hist = parent_hist.Clone()
    noise_subtraction_fit = TF1('noise_sub', Snoise, ped_upperlim, upper_fitlim, 4)
    noise_subtraction_fit.SetParameters(N3, Q0, alpha, mu_noise)

    #peaks_hist.Add(ped_subtraction_fit, -1)

    #peaks_hist.Add(noise_subtraction_fit, -1)

    # "N", "Q0", "Q1", "mu", "sigma1"
    peaks_ig = [800, Q0, 1.6, mu_ig, 0.2]
    peaks_hist, peaks_fit = Fit_SN(peaks_hist, ped_upperlim, noise_lowerlim, peaks_ig)

    # draw the peaks histogram and fit
    peaks_hist.Draw()
    input()

    Npeaks = peaks_fit.GetParameter(0)
    mu_peaks = peaks_fit.GetParameter(3)
    Q1 = peaks_fit.GetParameter(2)
    sigma1 = peaks_fit.GetParameter(4)
    Q0_peaks = peaks_fit.GetParameter(1)

    Npeaks_err = peaks_fit.GetParError(0)
    mu_peaks_err = peaks_fit.GetParError(3)
    Q1err = peaks_fit.GetParError(2)
    sigma1err = peaks_fit.GetParError(4)
    # ------------------------------------------------------------------------------------#

    # mu_peaks, sigma0, Q0, sigma1, Q1, Npeaks, alpha, W
    total_fit = TF1('total_fit', SParent, 0, upper_fitlim, 8)
    total_fit.SetParNames('mu', 'sigma0', 'Q0', 'sigma1', 'Q1', 'N', 'alpha', 'W')
    total_fit.SetParameters(4, 0.2, 1, 1, 5, 500, 0.1, 0.1)

    total_fit.SetParLimits(0, 0.001, 8)  # mu
    total_fit.SetParLimits(1, 0.001, 1)  # sigma0
    total_fit.SetParLimits(2, 0.001, 1)  # Q0
    total_fit.SetParLimits(3, 0.001, 3)  # sigma1
    total_fit.SetParLimits(4, 0.001, 3)  # Q1
    total_fit.SetParLimits(5, 100, n)  # Npeaks
    total_fit.SetParLimits(6, 0.001, 1)  # Alpha
    total_fit.SetParLimits(7, 0, 1)  # W
    # total_fit.SetParLimits(8, 1, n/2)
    '''
    mu = total_fit.GetParameter(0)
    sigma0 = total_fit.GetParameter(1)
    Q0 = total_fit.GetParameter(2)
    sigma1 = total_fit.GetParameter(3)
    Q1 = total_fit.GetParameter(4)
    N_peaks = total_fit.GetParameter(5)
    alpha = total_fit.GetParameter(6)
    W= total_fit.GetParameter(7)


    mu_err = total_fit.GetParError(0)
    sigma0_err = total_fit.GetParError(1)
    Q0_err = total_fit.GetParError(2)
    sigma1_err = total_fit.GetParError(3)
    Q1_err = total_fit.GetParError(4)
    N_peaks_err = total_fit.GetParError(5)
    alpha_err = total_fit.GetParError(6)
    W_err = total_fit.GetParError(7)
    '''
    print(Npeaks)
    print(N)
    print(N3)

    # --------------------------------------------------------#
    # Stotal(x, par):	#mu, sigma0, Q0, sigma1, Q1, N, alpha, N_ped, Nnoise
    # mu_peaks, sigma0, Q0, sigma1, Q1, Npeaks, alpha, W
    total_fit2 = TF1('total_fit2', Stotal, 0, upper_fitlim, 12)
    total_fit2.SetParNames('mu', 'sigma0', 'Q0', 'sigma1', 'Q1', 'N', 'alpha', 'N_ped', 'Nnoise', 'Q0noise', 'muNoise')

    total_fit2.SetParameters(mu_peaks, sigma0, Q0, sigma1, Q1, Npeaks, alpha, N, N3, Q0_noise, mu_noise)

    # draw the total fit

    parent_hist.SetFillColor(7)

    gStyle.SetFitFormat('5.3f')
    gStyle.SetOptFit(1110)
    gStyle.SetOptStat(11)
    # parent_hist.Fit('total_fit', 'S')
    parent_hist.Draw('e')
    total_fit2.Draw('same')

    xx = 0
    chi2temp = 0
    chi2 = 0

    for i in range(0, bins, 1):
        xx = parent_hist.GetBinCenter(i)
        chi2temp = ((parent_hist.GetBinContent(i) - total_fit2.Eval(xx))) - parent_hist.GetBinError(i)
        chi2 += chi2temp * chi2temp

    pt = TPaveText(8, 400, 14, 600)
    pt.AddText("chi2/ndf = {:.2f}".format(chi2 / (bins - 8)))
    pt.AddText("mu = {:.3f} #pm {:.3f}".format(mu_peaks, mu_peaks_err))
    pt.AddText("Q0 = {:.3f} #pm {:.3f}".format(Q0, Q0err))
    pt.AddText("#sigma0 = {:.3f} #pm {:.3f}".format(sigma0, sigma0err))
    pt.AddText("Q1 = {:.3f} #pm {:.3f}".format(Q1, Q1err))
    pt.AddText("#sigma1 = {:.3f} #pm {:.3f}".format(sigma1, sigma1err))
    pt.Draw()

    input()
    # counts, bin_width, ignore = plt.hist(data, bins=bins, range=(-1, 15))
    # bins = 0.5*(bin_width[1:] + bin_width[:-1])
    # data_list = [np.array(bins), np.array(counts)]
    # parameters = np.array([mu_peaks, sigma0, Q0, sigma1, Q1, Npeaks, alpha, N, N3])
    # plt.plot(bins, Stotal_python(data_list, parameters))
    #
    # plt.show()

    # str = f"../PMT Plots/PMT0/{pm}.pdf"
    # c1.Print(str)


if __name__ == '__main__':
    main()
