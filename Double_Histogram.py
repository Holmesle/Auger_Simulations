import matplotlib.pylab as plt
from matplotlib.pyplot import hist, viridis
import numpy as np


print('Start ----------------------------- ')

#v Variables
N = 1000
E1 = .371 + np.random.rand(N*2) * (240-.371)
E2 = .371 + np.random.rand(N) * (240-.371)


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
    plt.hist(norm1, weights = n1/d1, bins = 240, range = (0,240))
    plt.hist(norm2, weights = n2/d2, bins = nbins2, range = (0,240))
 
# create data
x = np.random.normal(size=50000)
y = x * 3 + np.random.normal(size=50000)
 
# Big bins
def Hist2D(E1,E2):
    ratio = len(E1)/len(E2)
    nbins2 = int(240*ratio)
    n1, b1 = np.histogram(E1, bins = 240,range = (0,240))
    n2, b2 = np.histogram(E2, bins = 240,range = (0,240))
    print(np.shape(E1))
    print(len(b1)-1,len(b2)-1)
    plt.hist2d(E1, E2, bins=(240, 240), cmap=plt.cm.jet)


fig1 = plt.figure(1)
Hist2D(x,y)
cb = plt.colorbar()
cb.set_label('Intensity')
# fig2 = plt.figure(2)
# DoubleHist(E1,E2)
plt.show()