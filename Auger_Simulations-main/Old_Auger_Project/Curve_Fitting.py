import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
import statistics as st
from scipy import signal

'''
Data
'''
A = np.random.normal(-6, 2, 1000)
B = np.random.normal(6, 2, 1000)

'''
Normalize
'''
n1, b1 = np.histogram(A, bins=100, range=(-10, 10))
n2, b2 = np.histogram(B, bins=100, range=(-10, 10))
n3, b3 = np.histogram(np.hstack([A, B]), bins=100, range=(-10, 10))
d1 = b1[1] - b1[0]
d2 = b2[1] - b2[0]
d3 = b3[1] - b3[0]

'''
Gaussian funcitons
'''

def guassian(x, a, u, s):
    return a*abs(x[1] - x[0])*np.exp(-((x-u)**2)/(2*s**2))/(np.sqrt(2*np.pi*s**2))

def double_guassian(x, a1, u1, s1, a2, u2, s2):
    return guassian(x, a1, u1, s1) + guassian(x, a2, u2, s2)


def Normal_Eqns(values, int0):
    a = values[0]
    u = values[1]
    s = values[2]
    peak = max(int0)
    mu = st.mean(int0)
    sigma = st.stdev(int0)
    print(peak,mu,sigma)
    return[a/(np.sqrt(2*np.pi*s**2)) - peak, u - mu,  s - sigma]

def normalize(A):
    n, b = np.histogram(A, bins=100, range=(-10, 10))
    d = b[1] - b[0]
    out = plt.hist((b[:-1] + b[1:])/2, weights=n /
                   d, bins=100, range=(-10, 10))
    return out


def single_guassian_curve(n,b):
    #initial guesss
    a, u, s = 1,0,1
    # root finding function
    sol = optimize.root(Normal_Eqns, [a,u,s], args=(
        n), method='lm')
    a,u,s = sol.x
    x = np.arange(-10,-10,0.01)
    print(a,u,s)
    y_curve = guassian(x, a, u, s)
    plt.plot(x,y_curve,color = 'black')


out = plt.hist(A,bins=100, range=(-10, 10))
single_guassian_curve(out[0],out[1])
plt.show()


# kai squared, p value - better with  which guassian
# make up pairs, probability of particular orbital filling A shell
# code packages for probabilities. 

