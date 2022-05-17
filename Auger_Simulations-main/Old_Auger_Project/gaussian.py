import matplotlib.pyplot as plt
from numpy import exp, loadtxt, pi, sqrt
from Detector_Curve_Fitting import detector
from lmfit import Model

data = loadtxt('model1d_gauss.dat')
b, n, nbins, max_range = detector(int(20000), 0.005e-3, 'PP',
                                  'Pt195m', 'Si', [70, 100, 130], 20000*5, 2)
x = (b[:-1] + b[1:])/2
y = n

def gaussian(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (sqrt(2*pi) * wid)) * exp(-(x-cen)**2 / (2*wid**2))


gmodel = Model(gaussian)
result = gmodel.fit(y, x=x, amp=5, cen=5, wid=1)

print(result.fit_report())

plt.plot(x, y, 'o')
plt.plot(x, result.init_fit, '--', label='initial fit')
plt.plot(x, result.best_fit, '-', label='best fit')
plt.legend()
plt.show()
