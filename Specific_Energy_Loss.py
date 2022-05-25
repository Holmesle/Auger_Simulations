from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import Isotope_Source as IS
import Materials as MT

'''
This module requires initial energy inputs from Isotope_Source and material 
properties form Materials module and outputs the final and deposited 
energies, and track lengths to the Particle_Path module.

Isotope_Source -> Specific_Energy_Loss -> Energy_Deposited -> Detector
     Materials ->
'''

#### Declared in the Detector Module ####
n_decays = int(2)
isotope = 'I125'
material = 'Si'
m, E = IS.JANIS_properties(n_decays, isotope)
material_values = MT.material_dict(material)  # material values
n_rows = len(E[0,:])

'''
Specific Energy Loss
'''
def SpecificEnergyLoss(energy_matrix, material_values, s):
    energy_matrix = np.array(energy_matrix)
    # material Values
    k = material_values[0]  # material dep value
    J = material_values[1]  # mean ionization potential
    Z = material_values[2]  # atomic number
    A = material_values[3]  # atomic weight
    p = material_values[4]  # density (g/cc)
    # set zeros to 0.0001 to avoid runtime error
    if len(np.shape((energy_matrix))) > 1:
        energy_matrix[energy_matrix == 0] = 0.001
    # compute integral for specific energy loss, equation 6a (Joy & Luo 1989)
    SEL = (-785*p*Z)/(A*energy_matrix)*np.log(1.166*(energy_matrix + k*J)/J)*s
    if len(np.shape((energy_matrix))) > 1:
        SEL[energy_matrix == 0.001] = 0
    return np.array(SEL)

'''
Duality: 

s max when SEL = Ef, or when 0 = SEL - Ef so when SEL is min (i.e conservation of energy)
'''

# SEL is a logarithmic function of Ei so rootfinding is best approach

def f(s0,material_values,Ei):
    k = material_values[0]  # material dep value
    J = material_values[1]  # mean ionization potential
    Z = material_values[2]  # atomic number
    A = material_values[3]  # atomic weight
    p = material_values[4]  # density (g/cc)
    return Ei - (785*p*Z)/(A*Ei)* np.log(1.166*(Ei + k*J)/J)*s0

def maximize_s(material_values,E,s0):
    sol = optimize.root(f,[s0],args = (material_values,E), jac=None,method='hybr')
    s = sol.x           # !!! there is probably a beter way to gues s0, I cannot think of one
    SEL = SpecificEnergyLoss(E, material_values, s)
    return s, SEL

def max_s(E,n_rows,n_decays):
    # simulated values
    Ei = E.reshape(n_rows*n_decays)
    Ei = Ei[Ei!=0]
    s0 = np.exp(Ei/1e3)-1
    s, SEL= maximize_s(material_values, Ei,s0)
    unique_s = np.unique(s)
    unique_E = np.unique(Ei)
    # E_map_s = dict(zip(unique_E,unique_s)) probably a better way using a dictionary
    s_matrix = np.zeros(np.shape(E))
    for i in range(unique_E.size):
        s_matrix[E == unique_E[i]] = unique_s[i]
    return s_matrix

'''
Check with graphs
'''
def Duality_Graphs(E,n_rows,isotope, material):
    # simulated values
    Ei = E.reshape(n_rows*n_decays)
    Ei = Ei[Ei!=0]
    s0 = np.exp(Ei/1e3)-1
    s, SEL= maximize_s(material_values, Ei,s0)
    # theoretical values
    Ei_theory = np.linspace(np.min(E[E!=0]),np.max(E),1000)
    s0_theory = np.exp(Ei_theory/1e3)
    s_theory, SEL_theory= maximize_s(material_values, Ei_theory,s0_theory)
    # units
    E_units = 'keV'
    E_cf = 1e-3
    s_units = 'um'
    s_cf = 1e-4             # recall 1s = 1 Angstrom
    #graph style
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams.update({'font.size': 13})
    # -SEL or ED vs Ei
    # plt.title('Conservation of Energy: Initial Energy = -Specific Energy Loss')
    # plt.plot(Ei*E_cf,-SEL*E_cf,'bo', label = 'Simulated')
    # plt.plot(Ei_theory*E_cf,-SEL_theory*E_cf,'--', color ='black', label = 'Theory')
    # plt.xlabel(f'Intial Auger Electron Energy ({E_units})')
    # plt.ylabel(f'Specific Energy Loss ({E_units})')
    # plt.legend()
    # s vs Ei
    plt.figure()
    plt.plot(Ei_theory*E_cf,s_theory*s_cf,'--', color ='black', label = 'Theory')
    plt.plot(Ei*E_cf,s*s_cf,'ro',label = 'Simulated')
    # plt.plot(Ei*E_cf, s0,'go')
    # plt.plot(Ei_theory*E_cf,s0_theory,'--',color ='green')
    plt.title(f'x{len(Ei)} {isotope} Auger Electrons in a {material} plate')
    plt.xlabel(f'Initial Energy ({E_units})')
    plt.ylabel(f'Maximal Distance Traveled ({s_units})')
    plt.legend()
    # RE vs Ei
    # plt.figure()
    # plt.title('Conservation of Energy: Remaining Energy = 0')
    # plt.plot(Ei*E_cf,(Ei + SEL)*E_cf,'mo',label = 'Simulated')
    # # plt.plot(Ei_theory*E_cf,(Ei_theory + SEL_theory)*E_cf,'o', color ='black', label = 'Theory')
    # plt.xlabel(f'Ei ({E_units})')
    # plt.ylabel(f'RE ({E_units})')
    # plt.legend()
    plt.tight_layout()

# Duality_Graphs(n_rows,isotope,material)

'''
SEL and RE at constant distance
# const path length such that Ei - min(SEL) = 0
'''  
def Energy_Graphs(E,n_rows,isotope,material):
    Ei = E.reshape(n_rows*n_decays)
    Ei = Ei[Ei!=0]
    s0 = np.exp(Ei/1e3)-1
    s_const = np.min(maximize_s(material_values, Ei, s0)[0])
    E = E.reshape(n_rows*n_decays)
    E = np.array(E[E!=0])
    SEL = SpecificEnergyLoss(E, material_values, s_const)
    SEL = np.array(SEL[SEL!=0])
    # theoretical values
    E_theory = np.linspace(np.min(E),np.max(E),1000)
    E_const = np.min(E)
    SEL_theory = SpecificEnergyLoss(E_theory,material_values,s_const).reshape(1000)
    # Graph Design
    fig = plt.figure()
    fig.set_size_inches(8, 4)
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serif"] = "Times New Roman"
    plt.rcParams.update({'font.size': 13})
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)
    con_fact = 1e-3
    E_units = 'keV'
    # graph Ei vs SEL, with const s
    ax1.plot(E_theory*con_fact,-SEL_theory*con_fact,'--',color = 'black',label = 'Theory')
    ax1.plot(E*con_fact,-SEL*con_fact,'ro',label = "Simulated")
    # graph E vs RE, with constant s
    ax2.plot(E_theory*con_fact,(E_theory + SEL_theory)*con_fact,'--',color = 'black',label = 'Theory')
    ax2.plot(E*con_fact,(E + SEL)*con_fact,'ro',label = "Simulated")
    # graph title
    ax1.set_xlabel('Initial Electron Energy (keV)')
    ax2.set_xlabel('Initial Electron Energy (keV)')
    ax1.set_ylabel('Specific Energy Loss (keV/Å)')
    ax2.set_ylabel('Final Energy (keV)')
    plt.suptitle(f'x{len(Ei)} {isotope} Auger electrons in a {material} plate of width {round(s_const,2)}Å')
    plt.legend()
    plt.tight_layout()

# Duality_Graphs(E,n_rows,isotope, material)
# Energy_Graphs(E,n_rows,isotope,material)
# plt.show()



