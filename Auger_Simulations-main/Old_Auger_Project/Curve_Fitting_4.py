from matplotlib.pyplot import errorbar
import numpy as np
import matplotlib.pylab as plt
import statistics as st
from numpy.linalg.linalg import matrix_rank
from scipy import optimize
from scipy import signal
from Detector_Curve_Fitting import detector

'''
Data
'''

def normalize(A, bins, a,b):
    # plt.figure()
    n1, b1 = np.histogram(A, bins=bins, range=(a, b))
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    out = plt.hist(x, weights=y, bins=100, range=(-10, 10))
    return out

def normalize1(n1, b1, nbins):
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    out = plt.hist(x, weights= y , bins=nbins, range=(0, 1))
    return out

'''
Gaussian funcitons
'''

def guassian(x, a, u, s):
    return a*abs(x[1] - x[0])*np.exp(-((x-u)**2)/(2*s**2))/(np.sqrt(2*np.pi*s**2))

def double_guassian(x, a1, u1, s1, a2, u2, s2):
    return guassian(x, a1, u1, s1) + guassian(x, a2, u2, s2)

'''
Curve Fitting
'''
def curvefit1(n, b, a, u, s):
    # popt = optimal value parameters
    # pcov = covariance of popt
    # print('a = ', a, 'u = ', u, 's = ', s)
    popt, pcov = optimize.curve_fit(guassian, (b[1:] + b[:-1])/2, n, sigma=np.sqrt(n+0.5), absolute_sigma=True, p0=[a, u, s], method='lm', maxfev=50000)
    energy = (b[1:] + b[:-1])/2
    chisq = ((guassian(energy, *popt) - n)**2/(n+0.5)).sum()/(len(energy) - 3 - 1)
    plt.plot(b, guassian(b, *popt), '--', color='red', label = 'Single Gauss')
    return chisq

def curvefit2(n, b, a1, u1, s1, a2, u2, s2):
    popt, pcov = optimize.curve_fit(
        double_guassian, (b[1:] + b[:-1])/2, n, sigma=np.sqrt(n+0.5), absolute_sigma=True, p0=[a1, u1, s1, a2, u2, s2], method='lm', maxfev=50000)
    energy = (b[1:] + b[:-1])/2
    chisq = ((double_guassian(energy, *popt) - n)**2/(n+0.5)).sum()/(len(energy) - 6 - 1)
    plt.plot(b, double_guassian(b, *popt),'k', label = 'Double Gauss')
    return chisq

'''
Null Hypothesis
'''

def second_derivative(x, f):  # centered difference approximation
    h = x[0] - x[1]
    numerator = f[2:] - 2*f[1:-1] + f[:-2]
    f_prime = numerator/(h**2)
    return f_prime


def first_derivative(x, f):  # forward difference approximation
    h = x[0] - x[1]
    numerator = f[2:] - f[:-2]
    f_prime = numerator/(2*h)
    return f_prime

def optimal_curve(n, b):
    # single Gaussian inital guesses for a,s,u
    u = 0
    u = sum(n)/len(n)
    # ua = st.mean(n)
    s = np.sqrt(sum((n[:] - u)**2)/len(n))
    # sa = st.stdev(n)
    A = max(n)
    a = (np.sqrt(2*np.pi)*A*abs(b[1] - b[0]))
    # single Gaussian inital guesses for a,s,u
    #(1) find slope so we can find peaks
    b_mid = (b[1:] + b[:-1])/2
    derv1 = first_derivative(b_mid, n)
    derv2 = second_derivative(b_mid, n)
    critical = (derv1 == np.zeros(len(derv1)))
    maxima = (derv2 < np.zeros(len(derv1)))
    minima = (derv2 > np.zeros(len(derv1)))
    max_n = n[1:-1].copy()
    min_n = n[1:-1].copy()
    max_n[~critical] == 0
    min_n[~critical] == 0
    max_n[~maxima] == 0
    min_n[~minima] == 0
    max_b = b_mid[1:-1].copy()
    min_b = b_mid[1:-1].copy()
    max_b[max_n ==0] = 0
    min_b[max_n == 0] = 0
    # maxima = (derv1.all() == 0) and (derv2.all < 0)
    # minima = (derv1.all() == np.zeros(len(derv1))) and (derv2.all > np.zeros(len(derv2)))
    # print(b_mid[minima])
    # nsub = n[(100<b) and (b<750)]
    # chisq2 = curvefit2(n, b, a1, u1, s1, a2, u2, s2)
    # chisq1 = curvefit1(n, b, a, u, s)
    chisq2 = 0
    chisq1 = 0
    return chisq1, chisq2,maxima,minima

'''
Test Data
'''

chisq1_array = []
chisq2_array = []
peak_distances = []

E = np.linspace(3,4,2)
for i in range(0,len(E)):
    pA = E[i]
    pB = -E[i]
    d = -E[i]
    plt.figure(figsize = (4, 4))
    A = np.random.normal(E[i], 2, 1000)
    B = np.random.normal(d, 2, 1000)
    out = plt.hist(np.hstack([A, B]), bins=100, range=(-10, 10))
    n = out[0]
    b = out[1]
    # smooth out data with a filter
    res = abs(pB - pA)
    print('res = ', res)
    filter = guassian(np.linspace(-res*10, res*10+1), 1, 0, res)
    smooth_n = signal.fftconvolve(n, filter, 'same')
    chisq1, chisq2,maxima,minima = optimal_curve(smooth_n, b)
    chisq1_array.append(chisq1)
    chisq2_array.append(chisq2)
    peak_distance = abs(d - E[i])
    peak_distances.append(peak_distance)
    b_mid = (b[1:] + b[:-1])/2
    # print(np.shape(minima), np.shape(maxima), np.shape(n[1:-2]),np.shape(b[1:-1]))
    # plt.plot(n[1:-2][minima], b[1:-1][minima],
    #          'o', color='blue', label='crit')
    # plt.plot(n[1:-2][maxima], b[1:-1][maxima],
    #          'o', color='red', label='max')
    plt.legend()
    plt.xlabel('Energy (keV)')
    plt.ylabel('Intensity')
    plt.title(f'Curve Fitting at peak distance {peak_distance}')
    
# plt.figure(figsize=(4,4))
# plt.title('Peak Distance vs Reduced Chi^2')
# plt.xlabel('Peak Distance (keV)')
# plt.ylabel('Reduced Chi^2')
# plt.plot(peak_distances, chisq1_array, 'o', color='red', label='Single Gauss')
# plt.plot(peak_distances, chisq2_array, 'o',
#          color='black', label='Double Gauss')
# plt.legend()
plt.show()

'''
Experimental Data
'''

# chisq1_ex_array = []
# chisq2_ex_array = []
# peak_ex_distances = []

# for res in range(70, 110, 10):
#     print('res = ', res)
#     plt.figure(figsize=(3, 3))
#     b, n, nbins, max_range = detector(
#         int(20000), 0.005e-3, 'PP', 'Pt195m', 'Si', res, 20000*5, 2)
#     optimal_curve(n, b)
#     out = plt.hist(n, bins=100, range=(0, 1))
#     chisq1, chisq2 = optimal_curve(n, b)
#     print('single chisq', chisq1)
#     print('double chisq', chisq2)
#     chisq1_ex_array.append(chisq1)
#     chisq2_ex_array.append(chisq2)
#     # peak_ex_distance = c*2
#     # peak_ex_distances.append(peak_ex_distance)
#     plt.legend()
#     plt.xlabel('Energy (keV)')
#     plt.ylabel('Intensity')
#     plt.title(f'Curve Fitting at resolution: {res}')

# # plt.figure(figsize=(4, 4))
# # plt.title('Peak Distance vs Reduced Chi^2')
# # plt.xlabel('Peak Distance')
# # plt.ylabel('Reduced Chai^2')
# # plt.plot(peak_distances, chisq1_array, 'o', color='red', label='Single Gauss')
# # plt.plot(peak_distances, chisq2_array, 'o',
# #          color='black', label='Double Gauss')
# # plt.legend()
# plt.show()


# for res in range(70,110,10):
#     print('res = ', res)
#     b, n, nbins, max_range = detector(int(20000), 0.005e-3,'PP','Pt195m','Si', res, 20000*5, 2)
#     out = normalize1(n, b, 100)
#     optimal_curve(n, b)
# plt.show()
