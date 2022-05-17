import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import Isotope_Source as IS

#### Declared in the Detector Module ####
n_decays = 50
m, E_Pt195m = IS.JANIS_properties(n_decays, 'Pt195m')
m, E_I125 = IS.JANIS_properties(n_decays, 'I125')
m, E_Sb119 = IS.JANIS_properties(n_decays, 'Sb119')
m, E_Hg197 = IS.JANIS_properties(n_decays, 'Hg197')

'''
Plate Dimensions
'''
height = 40e-3  # width (x) 40mm
width = 40e-3  # height (y)  40mm
source_displacement = 0.005e-3
diameter = 0.4e-3  # 0.4 mm   # detector diameter (z)
thickness = source_displacement*2  # thickness of the plate 0.02mm
distance_source_to_inner_plate2 = diameter + source_displacement
colormap = 'plasmas'

'''
Initailize Positions of Thin Films & Trajectory Angles
'''
# note: source is always located at z = 0
# randomize x,y,z,phi,theta matrices n = n_electrons, m = n_decays

def Initial_Conditions(E,height,width):
    # energy information
    N = len(E.flatten())   # number of electrons including 0s
    dim = np.shape(E)  # dim = [n_decays, n_electrons]
    E_keV = E.flatten()/1000   # Energies in keV including zeros
    tot_electrons = len(E_keV[E_keV != 0])  # total electrons excluding zeros
    # randomize initial positions and trajectory angles
    x0 = (np.random.rand(N) - (np.ones(N)*0.5))*(np.ones(N)*(width))
    y0 = (np.random.rand(N) - (np.ones(N)*0.5))*(np.ones(N)*(height))
    z0 = np.zeros(N)
    phi = (np.random.rand(N)) * np.pi         # test angle distribution
    theta = (np.random.rand(N))*2*np.pi
    # reshape to the same dim as Energy
    pos_matrices = np.vstack([np.reshape(x0, dim), np.reshape(y0, dim), np.reshape(z0, dim)])
    angle_matrices = np.vstack([np.reshape(phi, dim), np.reshape(theta, dim)])
    pos = np.vstack([x0,y0,z0])
    angle = np.vstack([phi,theta])
    return pos, angle, pos_matrices,angle_matrices


'''
2D plot of Initial Positions on Thin Film
'''
def xy_scatter(E,colormap):
    cm = mpl.cm.get_cmap(colormap)
    # energy information
    E_keV = E.flatten()/1000   # Energies in keV including zeros
    tot_electrons = len(E_keV[E_keV != 0])  # total electrons excluding zeros
    # initial conditions
    pos = Initial_Conditions(E, height, width)[0]
    # plots
    plt.figure()
    cm = mpl.cm.get_cmap('plasma')
    sc = plt.scatter(pos[0][E_keV != 0]*1e3, pos[1][E_keV != 0]* 1e3, c=E_keV[E_keV != 0], cmap=cm)
    cbar = plt.colorbar(sc)
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    cbar.set_label('Energy (keV)')
    plt.title(f'Electrons = {tot_electrons}')

'''
Polar Plot of Angles
'''

def polar_scatter(E,colormap):
    cm = mpl.cm.get_cmap(colormap)
    # energy information
    N = len(E.flatten()) 
    E_keV = E.flatten()/1000 
    tot_electrons = len(E_keV[E_keV != 0]) 
    # initial conditions
    angle = Initial_Conditions(E, height, width)[1]
    # plots
    tot_electrons = len(E_keV[E_keV != 0]) # total electrons excluding zeros
    fig = plt.figure()
    ax = fig.add_subplot(projection = 'polar')
    r = np.ones(N)
    c = ax.scatter(angle[0][E_keV != 0], r[E_keV != 0], c=E_keV[E_keV != 0], cmap=cm)
    cbar = plt.colorbar(c)
    plt.title(f'Electrons = {tot_electrons}')
    cbar.set_label('Energy (keV)')
    fig = plt.figure()
    ax = fig.add_subplot(projection='polar')
    r = np.ones(N)
    c = ax.scatter(angle[1][E_keV != 0], r[E_keV != 0], c=E_keV[E_keV != 0], cmap=cm)
    cbar = plt.colorbar(c)
    plt.title(f'Electrons = {tot_electrons}')
    cbar.set_label('Energy (keV)')

'''
Historgram x, y distributions
'''

def xy_hist(E,c):
    # energy information
    E_keV = E.flatten()/1000   # Energies in keV including zeros
    tot_electrons = len(E_keV[E_keV != 0])  # total electrons excluding zeros
    # initial conditions
    pos = Initial_Conditions(E, height, width)[0]
    # plots
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.hist(pos[0][E_keV != 0]*1e3,100, color=c)
    ax2.hist(pos[1][E_keV != 0]*1e3,100, color=c)
    ax1.set_xlabel('x (mm)')
    ax2.set_xlabel('y (mm)')
    ax1.set_ylabel('Frequency')
    plt.title(f'Electrons = {tot_electrons}')
    plt.tight_layout()

'''
Historgram theta, phi distributions
'''


def polar_hist(E,c):
    # energy information
    E_keV = E.flatten()/1000
    tot_electrons = len(E_keV[E_keV != 0])
    # initial conditions
    angle = Initial_Conditions(E, height, width)[1]
    # plots
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.hist(angle[0][E_keV != 0], 100, color = c)
    ax2.hist(angle[1][E_keV != 0], 100, color=c)
    ax1.set_xlabel('phi (rad)')
    ax2.set_xlabel('theta (rad)')
    ax1.set_ylabel('Frequency')
    plt.title(f'Electrons = {tot_electrons}')
    plt.tight_layout()


'''
Plots
'''

# xy_scatter(E_I125,'plasma')
# polar_scatter(E_I125,'plasma')
# xy_hist(E_I125,'blue')
# polar_hist(E_I125,'blue')
# plt.show()

'''
Return Values
'''



