from matplotlib.pyplot import errorbar
import numpy as np
import matplotlib.pylab as plt
import statistics as st
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
def optimal_curve(n,b):
    # single Gaussian inital guesses for a,s,u
    u = sum(n)/len(n)
    ua = st.mean(n)
    s = np.sqrt(sum((n[:] - u)**2)/len(n))
    sa = st.stdev(n)
    A = max(n)
    a = (np.sqrt(2*np.pi)*A*abs(b[1] - b[0]))
    # generate peaks at +/- 1 stdev from the mean
    # this method only works if the peaks are not too far from one another.
    u1 = u - s/2
    u2 = u + s/2
    print('u1 = ',u1)
    print('u2 = ',u2)
    print('s = ', u)
    print('sa = ', ua)
    # find a and s for double gaussian
    n1 = n[n <= u]
    n2 = n[n >= u]
    a1 = (np.sqrt(2*np.pi)*max(n1)*abs(b[1] - b[0]))
    a2 = (np.sqrt(2*np.pi)*max(n2)*abs(b[1] - b[0]))
    s1 = st.stdev(n1)
    s2 = st.stdev(n2)
    # nsub = n[(100<b) and (b<750)] 
    chisq2 = curvefit2(n, b, a1, u1, s1, a2, u2, s2)
    chisq1 = curvefit1(n, b, a, u, s)
    return chisq1, chisq2

'''
Test Data
'''

chisq1_array = []
chisq2_array = []
peak_distances = []

E = np.arange(2,5,0.5)
for i in range(0,len(E)):
    print('c = ',E[i])
    d = -E[i]
    plt.figure(figsize = (4, 4))
    A = np.random.normal(E[i], 2, 1000)
    B = np.random.normal(d, 2, 1000)
    out = plt.hist(np.hstack([A, B]), bins=100, range=(-10, 10))
    chisq1, chisq2 = optimal_curve(out[0], out[1])
    # print('single chisq', chisq1)
    # print('double chisq', chisq2)
    chisq1_array.append(chisq1)
    chisq2_array.append(chisq2)
    peak_distance = abs(d - E[i])
    peak_distances.append(peak_distance)
    plt.legend()
    plt.xlabel('Energy (keV)')
    plt.ylabel('Intensity')
    plt.title(f'Curve Fitting at peak distance {peak_distance}')
    
plt.figure(figsize=(4,4))
plt.title('Peak Distance vs Reduced Chi^2')
plt.xlabel('Peak Distance (keV)')
plt.ylabel('Reduced Chi^2')
plt.plot(peak_distances, chisq1_array, 'o', color='red', label='Single Gauss')
plt.plot(peak_distances, chisq2_array, 'o',
         color='black', label='Double Gauss')
plt.legend()
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
