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

def normalize(A, bins, a, b):
    # plt.figure()
    n1, b1 = np.histogram(A, bins=bins, range=(a, b))
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    # out = plt.hist(x, weights=y, bins=100, range=(-10, 10))
    return out


def normalize1(n1, b1, nbins):
    d1 = b1[1] - b1[0]
    x = (b1[:-1] + b1[1:])/2
    y = n1/d1
    # out = plt.hist(x, weights=y, bins=nbins, range=(0, 1))
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
    popt, pcov = optimize.curve_fit(guassian, (b[1:] + b[:-1])/2, n, sigma=np.sqrt(
        n+0.5), absolute_sigma=True, p0=[a, u, s], method='lm', maxfev=50000)
    energy = (b[1:] + b[:-1])/2
    print('n', n)
    chisq = ((guassian(energy, *popt) - n)**2 /
             (n+0.5)).sum()/(len(energy) - 3 - 1)
    plt.plot(b, guassian(b, *popt), '--', color='red', label='Single Gauss')
    return chisq


def curvefit2(n, b, a1, u1, s1, a2, u2, s2):
    popt, pcov = optimize.curve_fit(
        double_guassian, (b[1:] + b[:-1])/2, n, sigma=np.sqrt(n+0.5), absolute_sigma=True, p0=[a1, u1, s1, a2, u2, s2], method='lm', maxfev=50000)
    energy = (b[1:] + b[:-1])/2
    chisq = ((double_guassian(energy, *popt) - n) **
             2/(n+0.5)).sum()/(len(energy) - 6 - 1)
    plt.plot(b, double_guassian(b, *popt), '--',color = 'lime', label='Double Gauss')
    return chisq


'''
Null Hypothesis
'''

def optimal_curve(n, b, res):
    filter = guassian(np.linspace(-res*10, res*10+1), 1, 0, res)
    smooth_n = signal.fftconvolve(n, filter, 'same')
    b_mid = (b[1:] + b[:-1])/2
    # single Gaussian inital guesses for a,s,u
    u = 0
    u = sum(n)/len(n)
    # ua = st.mean(n)
    s = np.sqrt(sum((n[:] - u)**2)/len(n))
    # sa = st.stdev(n)
    A = max(n)
    a = (np.sqrt(2*np.pi)*A*abs(b[1] - b[0]))
    # double Gaussian inital guesses for a,s,u
    #(1) find slope so we can find peaks
    peaks, properties = signal.find_peaks(smooth_n, height = 20, distance = res)
    peak_n = smooth_n[peaks]
    peak_b = b_mid[peaks]
    smooth_n_non_zero = smooth_n[smooth_n!=0]
    begining = np.where(smooth_n == smooth_n_non_zero[0])
    end = np.where(smooth_n == smooth_n_non_zero[-1])
    b_begining = b[begining]
    b_end =b[end]
    s1 = (peak_b[0] - b_begining)[0]
    s2 = (b_end - peak_b[-1])[0]
    a1 = peak_n[0]*2*np.pi
    a2 = peak_n[-1]*2*np.pi
    u1 = peak_b[0]
    u2 = peak_b[-1]
    # print(a1, u1, s1, a2, u2, s2)
    # nsub = n[(100<b) and (b<750)]
    chisq2 = curvefit2(n, b, a1, u1, s1, a2, u2, s2)
    chisq1 = curvefit1(n, b, a, u, s)
    return chisq1, chisq2, peak_n, peak_b

'''
Test Data
'''

def test():  
    chisq1_array = []
    chisq2_array = []
    peak_distances = []

    E = np.linspace(2, 5, 5)
    for i in range(0, len(E)):
        plt.figure(figsize=(4, 4))
        pA = E[i]
        pB = -E[i]
        A = np.random.normal(pA, 2, 1000)
        B = np.random.normal(pB, 2, 1000)
        out = plt.hist(np.hstack([A, B]), bins=100, range=(-10, 10))
        n = out[0]
        b = out[1]
        # smooth out data with a filter
        res = abs(pB - pA)
        print('res = ',res)
        # filter = guassian(np.linspace(-res*10, res*10+1), 1, 0, res)
        # smooth_n = signal.fftconvolve(n, filter, 'same')
        # smooth_n = 1*n
        b_mid = (b[1:] + b[:-1])/2
        # plt.plot(b_mid,smooth_n,'--',color = 'black')
        chisq1, chisq2, peak_n, peak_b = optimal_curve(n, b, res)
        chisq1_array.append(chisq1)
        chisq2_array.append(chisq2)
        peak_distances.append(res)
        plt.plot(peak_b,peak_n,'or')
        plt.xlabel('Energy (keV)')
        plt.ylabel('Intensity')
        plt.xlim(-20,20)
        plt.ylim(0, max(peak_n)+10)
        plt.title(f'Curve Fitting at peak distance {res}')

    plt.figure(figsize=(4,4))
    plt.title('Peak Distance vs Reduced Chi^2')
    plt.xlabel('Peak Distance (keV)')
    plt.ylabel('Reduced Chi^2')
    plt.plot(peak_distances, chisq1_array, 'o', color='red', label='Single Gauss')
    plt.plot(peak_distances, chisq2_array, 'o', color='black', label='Double Gauss')
    plt.legend()
    plt.show()

'''
Experimental Data
'''

def experimental():
    chisq1_array = []
    chisq2_array = []
    peak_distances = []

    for res in range(80, 120,10):
        plt.figure(figsize=(4, 4))
        n, b = detector(int(20000), 0.005e-3, 'PP', 'Pt195m', 'Si', 100, 20000*5, 2)
        out = plt.hist(b, weights=n, bins=100, range = [0,1])
        n = out[0]
        b = out[1]
        print(np.shape(n), np.shape(b))
        # smooth out data with a filter
        print('res = ',res)
        filter = guassian(np.linspace(-res*10, res*10+1), 1, 0, res)
        smooth_n = signal.fftconvolve(n, filter, 'same')
        b_mid = (b[1:] + b[:-1])/2
        # plt.plot(b_mid,smooth_n,'--',color = 'black')
        chisq1, chisq2, peak_n, peak_b = optimal_curve(smooth_n, b, res)
        chisq1_array.append(chisq1)
        chisq2_array.append(chisq2)
        peak_distances.append(res)
        plt.plot(peak_b,peak_n,'or')
        plt.xlabel('Energy (keV)')
        plt.ylabel('Intensity')
        plt.ylim(0, max(peak_n)+10)
        plt.title(f'Curve Fitting at peak distance {res}')

    plt.figure(figsize=(4,4))
    plt.title('Peak Distance vs Reduced Chi^2')
    plt.xlabel('Peak Distance (keV)')
    plt.ylabel('Reduced Chi^2')
    plt.plot(peak_distances, chisq1_array, 'o', color='red', label='Single Gauss')
    plt.plot(peak_distances, chisq2_array, 'o', color='black', label='Double Gauss')
    plt.legend()
    plt.show()


# test()
experimental()
