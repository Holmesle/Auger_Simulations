import numpy as np
import matplotlib.pylab as plt
'''
Initialize Constants/Parameters
'''
hbar = 1    # hbar eV*a0
m = 1       # mass eV/c^2
k = 100
# the choice of units makes the energy unit the Hartree 27.2eV
# and the units of length the atomic bohr radius a0
p = hbar *k
Vmax = 0.005

'''
Initialize Dimensions
'''
L = 10
T = 100
dx = 0.001
eps = 0.5
dt = eps*(hbar**2/(2*m*dx**2 + Vmax))**-1 # to avoid numerical instability
x = np.arange(0,L,dx)
t = np.arange(0,T,dt)
s = 0.5

'''
Initialize Psi Matrices
'''
re_psi = np.zeros([len(x),len(t)])
im_psi = np.zeros([len(x),len(t)])

'''
Wave Packet Function
'''

def wavepacket(x,x0,s,k,complex):
    Gauss = np.exp(-(x-x0)**2/(2*s**2))*np.exp(-(x-x0)**2/(2*s**2))
    if complex == True:
        Plane = np.sin(k*x - x0)
    else:
        Plane = np.cos(k*x - x0)
    return Gauss*Plane

'''
Potential Function
'''
def V(x):
    return 0.5*x**2

'''
Initial Conditions
'''

x0 = 5
re_psi[:,0] = wavepacket(x,x0,s,k,False)
re_psi[:, 1] = wavepacket(x, (x0 + p/m*dt), s, k, False)
im_psi[:, 0] = wavepacket(x, x0, s, k, True)
im_psi[:, 1] = wavepacket(x, (x0 + p/m*dt), s, k, True)

'''
Inital Conditions Plots
'''

# fig = plt.figure(figsize=[10,6])
# ax1 = fig.add_subplot(2,1,1)
# ax2 = fig.add_subplot(2,1,2)
# ax1.plot(x,re_psi[:,0], color = 'blue', label = 'RE')
# ax2.plot(x,re_psi[:,1],color = 'blue', label = 'RE')
# ax1.plot(x, im_psi[:, 0], color='orange', label='IM')
# ax2.plot(x, im_psi[:, 1], color='orange', label='IM')
# ax1.legend()
# ax2.legend()
# plt.show()

'''
Future Wave Function Loops
'''
c1 = hbar*dt/(2*m*dt**2)
c2 = dt*hbar

for j in range(0,len(t)-1):
    for i in range(1,len(x)-1):
        re_psi[i, j+1] = re_psi[i, j]
        - c1*(im_psi[i+1, j] + im_psi[i - 1, j] - 2*im_psi[i, j+1])
        + c2*V(x)*re_psi[i, j+1]
    for i in range(1,len(x)-1):
        im_psi[i,j+1] = im_psi[i,j]
        + c1*(re_psi[i+2,j] + re_psi[i - 2,j] - 2*re_psi[i,j])
        - c2*V(x)*re_psi[i,j] 

print(im_psi[0:5,0:5])
print(re_psi[0:5, 0:5])





