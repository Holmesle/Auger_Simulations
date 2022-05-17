import numpy as np
import matplotlib.pyplot as plt
from numpy.core.records import array

'''
ELECTROSTATIC POTENTIAL
COMPUTATIONAL PHYSICS 366
WINTER 2022

This program ...


1. POTENTIAL FUNCTION

Computes the potential at any point in space, given as inputs,
and outputs the potential caused by four equal point charges
forming a square on an xy-plane. 
'''
# set length unit (L) for reference of distance between potentials
# where each point change is 2L distance from each other.
L = 0.001 # meter(s)
q = 1e-10
multipole_positions = np.array([np.array([L,L,-L,-L]), np.array([L,-L,-L,L]),np.array([0,0,0,0])])
multipole_magnitudes = np.array([q,q,q,q]) # Coloumb(s)
k = 8.99*10**9 # Coloumbs const (N*m^2/C^2)
N = 101  # number of values in the field point position arrays
# L_array = np.linspace(-2*L, 2*L, N) # field point positions
L_array = np.arange(-2*L, 2*L, 4*L/N+1)  # field point positions

def electrostatic_potential(x,y,z,x_prime,y_prime,z_prime,q):
    field_positions = [x,y,z]
    field_arrays = []
    position_types = np.array([type(x),type(y),type(z)])
    if type(x) or type(y) or type(z) == np.ndarray:
        # source values need to match shape(4,N)
        x_prime = np.transpose(np.array([x_prime]*N))
        y_prime = np.transpose(np.array([y_prime]*N))
        z_prime = np.transpose(np.array([z_prime]*N))
        q = np.transpose(np.array([q]*N))
        # field points need to match shape(4,N)
        # if not an array, the value need to be turned into array of length N
        for i in range(3):
            if position_types[i] != np.ndarray:
                # scalar values need to be N length arrays
                field_arrays.append(np.array([field_positions[i]*np.ones(N)]*4))
            else:
                # arrays need to be arrays of shape(4,N)
                field_arrays.append(np.array([field_positions[i]]*4))
        x,y,z = np.array(field_arrays)
    # compute difference vector (r) and then potential (V)
    r = np.sqrt((x - x_prime)**2 + (y - y_prime)**2 + (z - z_prime)**2)
    V = sum(k*q/r)
    return V


# print('V(x,y,z) = ', electrostatic_potential(
#     0, 0, 0, *multipole_positions, multipole_magnitudes), 'V')

'''
3. POWER SERIES APPROXIMATION

C
'''
def power(x,a):
    return((a)**-2)*(x - a)
    
'''
2. PLOT ELECTROSTATIC POTENTIAL

Create plots using the above function as a function of x, y, and z 
individually, while holding the other directions constant at 0. Each
plot extends an extra L distance past each point charge along the 
chosen direction.
'''

# # intiate figure
# fig = plt.figure()
# fig.set_size_inches(10, 5)
# ax1 = fig.add_subplot(1,3,1)
# ax2 = fig.add_subplot(1,3,2)
# ax3 = fig.add_subplot(1,3,3)

# conversion factors
m = 1e-3
c = 1000
unit1 = 'mm'
unit2 = 'kV'

# # X-axis
# ax1.plot(L_array/m, electrostatic_potential(
#     L_array, 0, 0, *multipole_positions, multipole_magnitudes)/c)
# ax1.set_xlabel(f'Position in x ({unit1})')
# ax1.set_ylabel(f'Electrostatic Potential ({unit2})')

# # Y-axis
# ax2.plot(L_array/m, electrostatic_potential(
#     0, L_array, 0, *multipole_positions, multipole_magnitudes)/c)  
# ax2.set_xlabel(f'Position in x ({unit1})')
# ax2.set_ylabel(f'Electrostatic Potential ({unit2})')

# # Z-axis
# ax3.plot(L_array/m, electrostatic_potential(
#     0, 0, L_array, *multipole_positions, multipole_magnitudes)/c)
# ax3.set_xlabel(f'Position in x ({unit1})')
# ax3.set_ylabel(f'Electrostatic Potential ({unit2})')

# plt.tight_layout()

# Heat map
fig1 = plt.figure()
fig1.set_size_inches(5, 5)
Z = V = electrostatic_potential(
    L_array, L_array, 0, *multipole_positions, multipole_magnitudes)/c
x = L_array
y = L_array
print(np.shape(x,y,Z))
plt.pcolormesh(x, y, Z)
plt.show()
y,x = np.mgrid[L_array,L_array]
c = plt.pcolormesh(L_array,L_array,V,cmap = 'plasma')
plt.colorbar
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()
plt.show()
