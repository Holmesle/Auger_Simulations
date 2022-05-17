import matplotlib.pylab as plt
from matplotlib.pyplot import hist, viridis
import numpy as np


print('Start ----------------------------- ')

#v Variables
N = 1000
# E1 = .371 + np.random.rand(N*2) * (240-.371)
# E2 = .371 + np.random.rand(N) * (240-.371)


# double histogram
def DoubleHist(E1,E2):
    ratio = len(E2)/len(E1)
    nbins2 = int(240*ratio)
    n1, b1 = np.histogram(E1, bins = 240,range = (0,240))
    n2, b2 = np.histogram(E2, bins = nbins2,range = (0,240))
    bins1 = (b1[1:] + b1[:-1])/2
    bins2 = (b2[1:] + b2[:-1])/2
    d1 = b1[1] - b1[0]
    d2 = b2[1] - b2[0]
    # plt.hist(bins1, weights = n1/d1, bins = 100, range = (0,10), label = 'x')
    # plt.hist(bins2, weights = n2/d2, bins = 100, range = (0,10), label = 'y')
    plt.hist(x, bins=100, range=(0, 10), label='x')
    plt.hist(y, bins=100, range=(0, 10), label='y')
    plt.xlabel('Bins')
    plt.ylabel('Intensity')
    plt.legend()
 
# create data
x = 5 + np.random.normal(size=50000)
y = 5 + np.random.normal(size=50000)*2
 
# Big bins
def Hist2D(E1,E2):
    ratio = len(E1)/len(E2)
    nbins2 = int(240*ratio)
    n1, b1 = np.histogram(E1, bins = 240,range = (0,240))
    n2, b2 = np.histogram(E2, bins = 240,range = (0,240))
    print(E1[:10],E2[:10])
    plt.hist2d(E1, E2, bins=(100, 100), range = ((0,10),(0,10)), cmap=plt.cm.jet)
    plt.scatter(E1,E2)
    plt.xlabel('x')
    plt.ylabel('y')


fig1 = plt.figure(1)
Hist2D(x,y)
cb = plt.colorbar()
cb.set_label('Intensity')
fig2 = plt.figure(2)
DoubleHist(x,y)
plt.show()
