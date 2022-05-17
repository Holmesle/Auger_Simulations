'''
This code computes the Spherical Pobabilities and plots them on 
a spherical plot. 
'''
print('---------------------------------')
print('NEW RUN')
print('---------------------------------')

import matplotlib.pylab as plt
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

# Varaibles

ℓ = 1
m = 0
n = 50
theta = np.linspace(0, 2*np.pi, n)
phi = np.linspace(0, np.pi, n)

print('phi=',len(phi))
print('theta=',len(theta))

# spherical to cartsian

def Sphere2Cart(phi,theta,r):
    x = r*np.sin(theta)*np.cos(phi,)
    y = r*np.sin(theta)*np.cos(phi,)
    z = r*np.cos(theta)
    return [x,y,z]


# 1. Spherical Harmonics

def SphericalHarmonics(ℓ, m, theta, phi):
    if ℓ == 0 and m == 0:
        return np.sqrt(1/(4*np.pi)) + 0*phi
    if ℓ == 1:
        if m == 1:
            return np.sqrt(3/(8*np.pi))*np.sin(theta)*(np.exp(1j*phi))
        if m == -1:
            return np.sqrt(3/(8*np.pi))*np.sin(theta)*(np.exp(-1j*phi))
        print(len(phi),len(theta))
        return np.sqrt(3/(4*np.pi))*np.cos(theta)
    if ℓ == 2:
        if m == 1:
            return np.sqrt(15/(8*np.pi))*np.sin(theta)*np.cos(theta)*(np.exp(m*1j*phi))
        if m == -1:
            return np.sqrt(15/(8*np.pi))*np.sin(theta)*np.cos(theta)*(np.exp(m*1j*phi))
        if m == 2:
            return np.sqrt(15/(32*np.pi))*np.sin(theta)**2*(np.exp(m*1j*phi))
        if m == -2:
            return np.sqrt(15/(32*np.pi))*np.sin(theta)**2*(np.exp(m*1j*phi))
        return np.sqrt(3/(8*np.pi))*np.sin(theta)*(np.exp(m*1j*phi))

# 2. Create a spherical plot function

def SpherePlot(ℓ, m, theta, phi):
    Y = SphericalHarmonics(ℓ, m, theta, phi)
    print(Y)
    print(len(phi),len(theta),len(Y))
    x, y, z = Sphere2Cart(phi, theta, Y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(x, y, z, linewidth = 0.5,cmap='plasma',edgecolor='none')
    plt.show()

SpherePlot(ℓ, m, theta, phi)
