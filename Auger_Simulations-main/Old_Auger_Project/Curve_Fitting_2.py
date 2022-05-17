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
A = np.random.normal(-5,1,1000)
B = np.random.normal(5,1,1000)
b,n, nbins, max_range = detector(int(20000), 0.005e-3, 'PP',
                'Pt195m', 'Si', [70, 100, 130], 20000*5, 2)

# def Range(nz,b):
#     b = (b[1:] + b[:-1])/2
#     # sort beginning and end points of peaks
#     As = b.copy()
#     Bs = b.copy()
#     print(np.shape(b),np.shape(nz))
#     filter1 = abs(nz[1:] - nz[:-1]) > np.ones(len(nz)-1)
#     print(np.shape(filter1))
#     filter2 = nz[1:] > nz[:-1]
#     np.shape(filter2)
#     print(print(np.shape(nz)))
#     As[filter1] = 0
#     Bs[filter1] = 0
#     As[~filter2] = 0
#     Bs[filter2] = 0
#     print(As, Bs)
#     return As,Bs

def normalize(A, bins, a,b):
    n1, b1 = np.histogram(A, bins=bins, range=(a, b))
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    R = x[n1 != 0]
    out = plt.hist(x, weights=y, bins=100, range=(R[0], R[-1]))
    return out

def normalize1(n1, b1, nbins,a,b):
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    R = x[n1 != 0]
    print(R[0], R[-1])
    print('nz = ',np.nonzero(n1))
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
    print('a = ', a, 'u = ', u, 's = ', s)
    popt, pcov = optimize.curve_fit(guassian, (b[1:] + b[:-1])/2, n, p0=[a, u, s], method='lm',maxfev = 2000)
    plt.plot(b, guassian(b, *popt),'--', color = 'red')
    # plt.title('single gauss')

def curvefit2(n, b, a1, u1, s1, a2, u2, s2):
    popt, pcov = optimize.curve_fit(
        double_guassian, (b[1:] + b[:-1])/2, n, p0=[a1, u1, s1, a2, u2, s2], method='lm')
    plt.plot(b, double_guassian(b, *popt),'k')
    # plt.title('double gauss')

'''
Null Hypothesis
'''
def optimal_curve(n,b):
    N = int(np.shape(n)[0]/2)
    print(max(n[:N]), max(n[N:]))
    # nsub = n[(100<b) and (b<750)]
    curvefit2(n, b, max(n[:N]), st.mean(n[:N]), st.stdev(
        n[:N]), max(n[N:]), st.mean(n[N:]), st.stdev(n[N:]))
    # curvefit1(n, b, max(n), st.mean(n), st.stdev(n))

def reduced_chi(out1):
    pass


# plt.figure()
# outs = plt.hist(np.hstack([A,B]), bins= 100, range=(-10, 10))
# nz =  list(np.nonzero(np.hstack([A, B])))
# nz = np.reshape(nz,np.shape(nz)[1])
# # print(np.shape(nz))
# # plt.show()
# Ranges = Range(nz, outs[1])
# # out = normalize1(np.hstack([A, B]), 100, Range[0][0], Range[1][0])
# # optimal_curve(out[0], out[1])

plt.figure()
# Ranges = Range(n[0],b[0])
out = normalize1(n[0], b[0], 100,-10,10)
optimal_curve(n[0], b[0])
plt.show()

# plt.figure()
# out = normalize1(n[1], b[1], 100)
# optimal_curve(n[1], b[1])
# plt.show()

# plt.figure()
# out = normalize1(n[2], b[2], 100)
# optimal_curve(n[2], b[2])
# plt.show()
