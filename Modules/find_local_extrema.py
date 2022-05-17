import numpy as np
import matplotlib.pyplot as plt

x_range = np.arange(-10,10,1)

'''
Derivatives - Central Difference
'''

def second_derivative(x, f):  # centered difference approximation
    return (f[2:] - 2*f[1:-1] + f[:-2])/((x[0] - x[1])**2)

def first_derivative(x, f):  # forward difference approximation
    return (f[2:] - f[:-2])/(2*(x[0] - x[1]))

'''
Local Extrema
'''
def Local_Extrema(x,f):
    # make coppies of arrays
    x_mid = (x_range[2:]+x_range[:-2])/2
    maxima, minima = x_mid.copy(), x_mid.copy()
    # find 2nd derivative = 0, extrema
    extrema =  second_derivative(x, f) == 0
    maxima[~extrema] = 0
    minima[~extrema] = 0
    # find 1st derivative +, i.e. increasing (1) or decreasing (0)
    positive_slope = first_derivative(x, f) > 0
    maxima[positive_slope] = 0
    minima[~positive_slope] = 0
    # return sorted list [value, index, min/max] from high to low
    index_max = np.nonzero(maxima)
    index_min = np.nonzero(minima)
    max_values, min_values = maxima[maxima != 0], minima[minima != 0]
    return [[max_values[:], index_max[:]], [min_values[:], index_min[:]]]
    
'''
Plots
'''

x_mid = (x_range[2:]+x_range[:-2])/2
y = x_range**2
y_prime = first_derivative(x_range, y)
y_double_prime = second_derivative(x_range, y)

print('max values:',Local_Extrema(x_range, y)[0][0])
print('min values:', Local_Extrema(x_range, y)[1][0])
print(x_range[10])

plt.plot(x_range,y,'ro', label = 'f(x)')
# plt.plot(x_mid, y_prime, label="f'(x)")
# plt.plot(x_mid, y_double_prime, label="f'(x)")
# plt.legend()
plt.show()
