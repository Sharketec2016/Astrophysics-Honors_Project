import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
from scipy.optimize import minimize

# def f0(x, sigma0, N, Q0, W):
#     return N*( (1-W)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)


# def f1(x, sigma1, N1, Q1, Q0):
#     return N1*( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2))


# def f2(x, sigma1, N1, Q1, Q0):
#     return N1*( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 

# def f3(x, sigma1, N1, Q1, Q0):
#     return N1*( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2))

# def f4(x, sigma1, N1, Q1, Q0):
#     return N1*( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2))

# def f5(x, sigma1, N1, Q1, Q0):
#     return N1*(((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))

# def f6(x, sigma1, N1, Q1, Q0):
#     return N1*(((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2))

# def S_noise(x, Q0, alpha, W, N):
#     return np.where(x<Q0, 0, N*alpha*W*np.exp(-alpha*(x-Q0) - mu))




def g(x, mu, sigma0, Q0, sigma1, Q1, N, W, alpha): 
    
    return np.where(x<Q0, N*(  ( (1-W)/(np.sqrt(2*np.pi)*sigma0) ) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)
            + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))
            + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1-Q0)/sigma1)**2))  
            + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1-Q0)/sigma1)**2)) ), 
            
            
            N*(  ( (1-W)/(np.sqrt(2*np.pi)*sigma0)) * np.exp(-0.5*((x-Q0)**2/(sigma0)**2) - mu)  + (alpha*W*np.exp(-alpha*(x-Q0) - mu))
            + ( ((mu**1 * np.exp(-mu))/sp.special.factorial(1)) * (1/(np.sqrt(2*np.pi)*1*sigma1)) * np.exp((-1/(2*1)) * ((x-1*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**2 * np.exp(-mu))/sp.special.factorial(2)) * (1/(np.sqrt(2*np.pi)*2*sigma1)) * np.exp((-1/(2*2)) * ((x-2*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**3 * np.exp(-mu))/sp.special.factorial(3)) * (1/(np.sqrt(2*np.pi)*3*sigma1)) * np.exp((-1/(2*3)) * ((x-3*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**4 * np.exp(-mu))/sp.special.factorial(4)) * (1/(np.sqrt(2*np.pi)*4*sigma1)) * np.exp((-1/(2*4)) * ((x-4*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**5 * np.exp(-mu))/sp.special.factorial(5)) * (1/(np.sqrt(2*np.pi)*5*sigma1)) * np.exp((-1/(2*5)) * ((x-5*Q1-Q0)/sigma1)**2))
            + ( ((mu**6 * np.exp(-mu))/sp.special.factorial(6)) * (1/(np.sqrt(2*np.pi)*6*sigma1)) * np.exp((-1/(2*6)) * ((x-6*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**7 * np.exp(-mu))/sp.special.factorial(7)) * (1/(np.sqrt(2*np.pi)*7*sigma1)) * np.exp((-1/(2*7)) * ((x-7*Q1-Q0)/sigma1)**2))  
            + ( ((mu**8 * np.exp(-mu))/sp.special.factorial(8)) * (1/(np.sqrt(2*np.pi)*8*sigma1)) * np.exp((-1/(2*8)) * ((x-8*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**9 * np.exp(-mu))/sp.special.factorial(9)) * (1/(np.sqrt(2*np.pi)*9*sigma1)) * np.exp((-1/(2*9)) * ((x-9*Q1-Q0)/sigma1)**2)) 
            + ( ((mu**10 * np.exp(-mu))/sp.special.factorial(10)) * (1/(np.sqrt(2*np.pi)*10*sigma1)) * np.exp((-1/(2*10)) * ((x-10*Q1-Q0)/sigma1)**2))) )    
    
    
    

    
    
    
    

# def error_bins(raw_data, fitted_data, raw_error):
#     chi2_table = {0: [0.00, 1.29], 1: [0.27, 2.75], 2: [0.74, 4.25], 3: [1.10, 5.30], 4: [2.34, 6.78], 5: [2.75, 7.81],
#                 6: [3.82, 9.28],
#                 7: [4.25, 10.30], 8: [5.30, 11.32], 9: [6.33, 12.79], 10: [6.78, 13.81], 11: [7.81, 14.82],
#                 12: [8.83, 16.29], 13: [9.28, 17.30],
#                 14: [10.30, 18.32], 15: [11.32, 19.32], 16: [12.33, 20.80], 17: [12.79, 21.81], 18: [13.81, 22.82],
#                 19: [14.82, 23.82], 20: [15.83, 25.30]}

#     for i in range(len(raw_data)):
#         if raw_data[i] < 10:
#             error_band = chi2_table.get(raw_data[i])
#             upper = abs(error_band[0] - raw_data[i])
#             lower = abs(error_band[1] - raw_data[i])
#             residual = fitted_data[i] - raw_data[i]
#             if(residual >= 0):
#                 raw_error[i] = upper
#             elif(residual < 0):
#                 raw_error[i] = lower
#     return raw_error
 
 
counts = np.array([  2,   5,   3,   3,  19,  45, 339, 307, 173, 240, 307, 436, 458, 396,
353, 387, 437, 423, 383, 356, 342, 338, 276, 295, 244, 233, 207, 187,
162, 165, 136, 125,  98,  97,  77,  93,  70,  64,  62,  58,  36,  44,
36,  35,  32,  36,  30,  35,  25,  29,  31,  26,  18,  29,  21,  28,
21,  19,  25,  21,  21,  18,  15,  31,  20,  25,  20,  20,  13,  26,])





bins = np.array([ 0.07857143,  0.23571429,  0.39285714,  0.55,        0.70714286,  0.86428571,
1.02142857,  1.17857143,  1.33571429,  1.49285714,  1.65      ,  1.80714286,
1.96428571,  2.12142857,  2.27857143,  2.43571429,  2.59285714,  2.75,
2.90714286,  3.06428571,  3.22142857,  3.37857143,  3.53571429,  3.69285714,
3.85      ,  4.00714286,  4.16428571,  4.32142857,  4.47857143,  4.63571429,
4.79285714,  4.95      ,  5.10714286,  5.26428571,  5.42142857,  5.57857143,
5.73571429,  5.89285714,  6.05      ,  6.20714286,  6.36428571,  6.52142857,
6.67857143,  6.83571429,  6.99285714,  7.15      ,  7.30714286,  7.46428571,
7.62142857,  7.77857143,  7.93571429,  8.09285714,  8.25      ,  8.40714286,
8.56428571,  8.72142857,  8.87857143,  9.03571429,  9.19285714 , 9.35,
9.50714286,  9.66428571,  9.82142857, 9.97857143, 10.13571429, 10.29285714,
10.45,       10.60714286, 10.76428571, 10.92142857])
 
def chisqfunc(arguments):
    chisq = np.sum(((counts - g(bins, *arguments))**2/counts))
    return chisq



new_counts=[]
new_bins=[]
for i in range(len(bins)):
    if(bins[i] < 7):
        new_bins.append(bins[i])
        new_counts.append(counts[i])

counts = new_counts
bins = new_bins



gain = 1.1
std = 0.25
gain0 = 1
std0 = .5
alpha = 0.001
W = 0.45
N = 2500
mu = 3.5
k = 1
#mu, sigma0, Q0, sigma1, Q1, N, W, alpha
initial_guess = [mu, std0, gain0, std, gain, N, W, alpha]


fig, ax = plt.subplots(figsize=(15, 10))



#here we are calculating the error on each bin. Since you need to have the model already done for it to calculate the error on it properally, then we use curve_fit to find optimal parameters to fill the model to then use for the chi2 minimization
bin_error = np.sqrt(counts)
ax.bar(bins, counts, edgecolor='black', width=.15, yerr=bin_error)



# cons = ({'type':'eq', 'fun':lambda x: x[8]%2}, {'type':'eq', 'fun':lambda x: x[8]%2+1})

#here we are using the curve_fit method that implements a least_squares method for minimzation to aquire some inital guesses for the chi2 minimization
popt, pcov = sp.optimize.curve_fit(g, bins, counts, p0=initial_guess, bounds=([0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001], [10, 10, 10, 10, 10, 7000, 10, 10]), sigma=bin_error)
plt.plot(bins, g(bins, *popt), label='Least Squares minimization', color='black')





#mu, sigma0, Q0, sigma1, Q1, N, W, alpha
bnds = ((0.001, 10), (0.001, 10), (0.001, 10), (0.001, 10), (0.001, 10), (0.001, 7000), (0.001, 10), (0.001, 10))#, (1, 20))
result =  minimize(chisqfunc, x0=popt, bounds=bnds)#, constraints=cons, method='trust-constr')
print(result)
for i in range(len(result.x)):
    print(result.x[i])

ax.plot(bins, g(bins, *result.x), color='orange', label='$\chi2$Minimication')






textstr = r'''$\mu$ = {:3f} +/- {:.3f}
$\sigma_0$ = {:.3f} +/- {:.3f}
$Q_0$ = {:.3f} +/- {:.3f}
$\sigma_1$ = {:.3f} +/- {:.3f}
$Q_1$ = {:.3f} +/- {:.3f}
N = {:.0f} +/- {:.3f}
W = {:.3f} +/- {:.3f}
$\alpha$ = {:.3f} +/- {:.3f}
$\chi2$/dof = ${:.1f}/63$'''.format(result.x[0],result.jac[0], result.x[1],result.jac[1], result.x[2],result.jac[2], result.x[3],result.jac[3],result.x[4],result.jac[4], result.x[5],result.jac[5], result.x[6],result.jac[6], result.x[7],result.jac[7], result.fun, 63)




# ax.plot(bins, f0(bins, result.x[1], result.x[5], result.x[2], result.x[6]), label=r'$S_{ped}$')
# ax.plot(bins, f1(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{1}$')
# ax.plot(bins, f2(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{2}$', color='black')
# ax.plot(bins, f3(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{3}$')
# ax.plot(bins, f4(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{4}$')
# ax.plot(bins, f5(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{5}$')
# ax.plot(bins, f6(bins, result.x[3], result.x[5], result.x[4], result.x[2]), label=r'$S_{6}$')

# # Q0, alpha, W, N
# ax.plot(bins, S_noise(bins, result.x[2], result.x[7], result.x[6], result.x[5]), label='S_noise')



props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(10-4, np.max(counts), textstr, fontsize=14,verticalalignment='top', bbox=props)
ax.set_xlabel(r'$[nC]$')
ax.set_ylabel('Histogram Counts')
ax.legend(fontsize=14)
plt.xticks([r for r in range(11)])
ax.grid()
plt.show()

