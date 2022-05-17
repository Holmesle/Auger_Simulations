import numpy as np
import matplotlib.pyplot as plt


# initialize variables
dx = 0.1
t_max = 50
x_max = 50

# calculate velocity from conditions
T = 1
mu = 2
v = np.sqrt(T/mu)

dt = 0.9*dx/v

# check numerical error

if dt >= dx/v:
    print('error: dx or v are not big enough')
print(v)
print(dx/v)
print(T/dt)

# x-position and time arrays
x = np.arange(0,x_max,dx)
t = np.arange(0,t_max,dt)

# y-position matrix
y = np.zeros([len(t), len(x)])

# Wave Packet function
def Gauss(x0,x,s):
    return np.sqrt(2*np.pi*s)*np.exp(-((x - x0)/2*s)**2)

# Initial Conditions
y[0,:] = Gauss(25,x,0.5)
y[1,:] = Gauss(25, x-v*dt, 0.5)

# Cause and Effect Loop
for i in range(1,len(t)-1):
    for j in range(1,len(x)-1):
        y[i+1,j] = 2*y[i,j] - y[i-1,j] 
        - (v*dt/dx)**2*(2*y[i,j]-y[i,j-1]-y[i,j+1])
    # plt.cla()  # Clear off what was previously drawn
    # plt.plot(x, y[i+1,:])
    # plt.xlim([0,50])
    # plt.ylim([-10,10])
    # plt.pause(0.1)  # Draw to the screen (and wait a moment)

# print(y[0,:])
# print(y[2,:] - y[1,:])

# animation
for i in range(len(t)):
    plt.cla()  # Clear off what was previously drawn
    plt.plot(x, y[i,:])
    plt.xlim([0,50])
    plt.ylim([-10,10])
    plt.pause(0.01)  # Draw to the screen (and wait a moment)
plt.show()

 