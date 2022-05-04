import numpy as np 
import scipy.fft
from scipy.fft import fft, fftfreq, fftshift
import cmath
import math

#Code to obtain the bias by performing the FFT setting the zero frequency to zero
#setting 95% of the highest frequencies to zero, then constructiong a linear relation for the integration of the points above zero
#subtracting the linear relation from the integration and finding where it goes to zero
#finally subtracting the point where the integration goes ti zero form the amplitude of the zero frequency of the FFT
#that value is the bias which is finally subtracted from the waveform.

def fft(volts):
    #function for78 doing the FFT setting the 95% of the highes freuqncies 
    # and the zero frequency to zero. Returning the inverse FFT. 

    n=len(volts) 
    Fft=scipy.fft.rfft(volts,n=n)
    #doing the FFT

    A=[]
    thetas=[]
    for j in Fft:
        r, theta=cmath.polar(j)
        A.append(r)
        thetas.append(theta)
    # for loop to go from complex nymber (X,Y) to its amplitude and angle regresentation (A,pi)

    entries=Fft.size
    percentage= math.floor(entries-(entries*0.95))
    #get the value where the 95% of the highest frequencies begin depending of the size of the FFT array

    for h in range (percentage, len(A)):
        A[h]=0
    #set the 95% of the highest frequencies to zero

        amplitude=A[0]/n
    # the amplitude of the zero frequecy to use is later

    A[0]=0
    #set amplitude of zero frequency to zero

    complex_n=[]
    for k in range (Fft.size):
        comp=cmath.rect(A[k],thetas[k])
        complex_n.append(comp)
    #go from (A,pi) to (X,Y) complex numbers

    rverFft=scipy.fft.irfft(complex_n, n=n)
    #reverse FFT

    return amplitude, rverFft
    #return the amplitude of the zero frequency and the reverse FFT


def bias (voltage,slope):
    #This function is to obtain the bias with the help pf the fft 
    #it finds the bias with a binary search, integration and linear fit
    #the inputs of this function are a list of voltages ( the waveform points) and the slope for the linear fit
    #which is the number of points in the waveform

    amplitude, reversFft=fft(voltage)
    #call the fft function

    ifft=np.array(reversFft)
    ifft_50=ifft+amplitude+50
    ifft_100=ifft+amplitude+100
    #list to calculate the two integration points for the linear fit
    ifft_50_integration=ifft_50.sum()
    ifft_100_integration=ifft_100.sum()
    #the two integration points to calculate the intercept for the linear relation

    b=((ifft_50_integration+ifft_100_integration)-(slope*(50+amplitude+100+amplitude)))/2
    #calculation of the intercept

    high=0
    low=50
    mid=0
    #highest point is zero because the y component of the 0 is the highest point
    #low component is a big number because the y component of that x components is really low
    #binary search to find the bias
    while low>=high:
        mid=(low+high)/2
        #get the average 
        # print(mid)
        increment=ifft+mid
        #add the average to the ifft so then we can integrate at that point by adding all the points
        #that are above zero
        y=(slope*mid)+b
        #get the linear with the parameters previouslly value related to the position of the integration point

        integration=0
        for e in increment:
            if (e>0):
                integration+=e
        #add every point that is above zero
        target=integration-y 
        #subtract the linear value to its respective integration point

        #binary search for a point within a range
        #the range choosen is from 10^-1 to 10^-2
        if (10**-2<=target and target<=10**-1):
            break
        #if the point related to the average is greater than 10^-2
        #and smaller than 10^-1 the point is found and returns the x component
        elif (target>10**-1):
            high=mid
            # if the point is greater that 10**-1 the high component equals to mid, in order to reduce options
        elif (target<10**-2):
            low=mid
        # if the point is grater that 10**-1 the high component equals to mid, in order to reduce options

    base=amplitude-mid
    base=amplitude
    # calculation of the baseline 
    
    return base

def voltage_base_and_max (volts, bias, time):
    #this function substracts the bias from the waveform
    #also returns the maximum amplitude and the time when this happens 
    #the inputs are the voltage array, its respective bias (calculated with the above functions) and the time spacing
 #   voltage=np.array(volts)

    #here is where the subraction is happening. Since i dont want to subract the negative parts off, b/c part of our waveform is negative, im negating the subtraction all together.
    bias=0
    voltage_base=volts-bias
    volts_base=voltage_base.tolist()

    volts_max=voltage_base.max()
    time_max_vb=time[volts_base.index(volts_max)]
    return voltage_base, volts_max, time_max_vb

