import numpy as np
import matplotlib.pyplot as plt

'''
Computational Integrating Charge Distribution
PH366 - Computational Physics
01/07/22
'''
'''
(0) Initiate Variables
'''
# ------------------------------------------- #
# Square Surface Charge
# ------------------------------------------- #
L = 1 # unit length (m)

# square of charge (square of length 2*L centered @ the origin)
dLp = L/50   # step size, for integration
Lp_array = np.arange(-L,L,dLp)
zp = np.zeros(len(Lp_array)) # plane of charge
sigma = 1  # surface charge density, C/m
k = 8.988*10**9  # kg*m^3*/(s^2*c^2)
k = 1 # kg*m^3*/(s^2*c^2)

# ------------------------------------------- #
# Field Points
# ------------------------------------------- #

dL = L/100    # step size, for resolution
L_array = np.arange(-2*L,2*L,dL) # one L past the square on each side

'''
(1) Potential function for at any point in space
'''
# takes a single point

# def V(x,y,z):
#     # assuming we are integrating via squares in the xy-plane (dx'dy')
#     # both are equal length arrays so we can use 1 for loop
#     r = 0
#     for dxp in Lp_array: #dx',dy'
#         for dyp in Lp_array:
#             r += 1/np.sqrt((x - dxp)** 2 + (y - dyp)**2 + (z - 0)**2)
#     V = (k*sigma)*dxp*dyp/r
#     return V


def V(x, y, z):
    I = 0
    dxp = Lp_array[1] - Lp_array[0]
    dyp = Lp_array[1] - Lp_array[0]
    for i in range(len(Lp_array)):  # dx',dy'
        for j in range(len(Lp_array)):
            I += (k*sigma)*dxp*dyp/np.sqrt((x - Lp_array[i])** 2 + (y - Lp_array[j])**2 + (z - 0)**2)
    V = (k*sigma)*dxp*dyp*I
    return V

# ------------------------------------------- #
# Test Function
# ------------------------------------------- #


xt, yt, zt = 0, 0, -1
Vt = V(xt, yt, zt)

v = 1  # unit convertsion for potential (V)
m = 1e-3  # unit conversion for length (mm)
c = 1  # unit conversion for charge (C)

print('\nThe potential of a uniform square charge distribution of density %1.2e C and width 2L (where L = %.2f mm)\
    \n at (x,y,z) = (%.2f mm,%.2f mm,%.2f mm) is %.2f V.' %(sigma/c,L/m,xt/m,yt/m,zt/m,Vt/v))

'''
(2) Plot the V(r) in each Dimension
'''
v = 1  # unit convertsion for potential
m = 1  # unit conversion for length
c = 1  # unit conversion for charge
m_unit = 'm'
v_unit = 'V'
c_unit = 'C'

fig = plt.figure()
fig.set_size_inches(12, 4)
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

# ------------------------------------------- #
# X Direction
# ------------------------------------------- #

ax1.plot(L_array/m, V(L_array, 0, 0)/v)
ax1.set_xlabel('Position in X ('+m_unit+')')
ax1.set_ylabel('Electrostatic Potential ('+v_unit+')')

# ------------------------------------------- #
# Y Direction
# ------------------------------------------- #

ax2.plot(L_array/m, V(0, L_array, 0)/v)
ax2.set_xlabel('Position in Y ('+m_unit+')')
ax2.set_ylabel('Electrostatic Potential ('+v_unit+')')

# ------------------------------------------- #
# Z Direction
# ------------------------------------------- #

ax3.plot(L_array/m, V(0, 0, L_array)/v)
ax3.set_xlabel('Position in Z ('+m_unit+')')
ax3.set_ylabel('Electrostatic Potential ('+v_unit+')')
ax3.set_xlim([-2*L-0.5,2*L+0.5])

plt.tight_layout()
plt.show()
