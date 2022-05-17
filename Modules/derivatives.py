import numpy as np
import matplotlib.pyplot as plt

x_range = np.arange(-2*np.pi,2*np.pi,0.01)

# def find_local_extrema(x, f):
#     h = x[1:] - x[:-1]
#     numerator = f[2:] - 2*f[] + f[]
#     return numerator/(h**2)

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


# def first_derivative(x, f):  # forward difference approximation
#     h = x[1:] - x[:-1]
#     numerator = f[1:] - f[:-1]
#     f_prime = numerator/(h*2)
#     return f_prime


x_mid = (x_range[2:]+x_range[:-2])/2
print(x_mid[:10],x_mid[-10:])
y = np.cos(x_range)
y_prime = first_derivative(x_range, y)
y_double_prime = second_derivative(x_range, y)
print(len(x_mid), len(y_prime), len(y))

plt.plot(x_range,y, label = 'f(x)')
plt.plot(x_mid, y_prime, label= "f'(x)")
plt.plot(x_mid, y_double_prime, label= "f'(x)")
plt.legend()
plt.show()
