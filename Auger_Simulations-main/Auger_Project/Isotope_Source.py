import matplotlib.pylab as plt
import numpy as np
import Janis_Data as Janis
import BRICC

'''
This module requires initial energy inputs from the BRICC and Janis Modules 
module and outputs the isotope source properites (mass, initial_energy_matrix) 
to the Main_Detector module.

     BRICC -> Isotope_Source -> Detector
Janis_Data -> 

'''


def JANIS_properties(n_decays,isotope):
    # mass of an electron in MeV/c^2 (wiki: https://en.wikipedia.org/wiki/Electron_rest_mass)
    mass = 0.511
    # generate energy matrix (shape = [n_decays,n_electrons], zeros are placeholders)
    initial_energy_matrix_eV = Janis.initial_energies(isotope, n_decays)
    return mass, initial_energy_matrix_eV


def BRICC_propterties(n_decays, isotope):
    # mass of an electron in MeV/c^2 (wiki: https://en.wikipedia.org/wiki/Electron_rest_mass)
    mass = 0.511
    # generate energy matrix (shape = [n_decays,n_electrons], zeros are placeholders)
    initial_energy_matrix_eV = BRICC.initial_energies(isotope, n_decays)
    return mass, initial_energy_matrix_eV

'''
Checks/Plots
'''
# Test parameters (declared in detector module)

decay_per_interval = 100000
n_intervals = 2
# Energy matrix = [n_decays,n_intervals]
n_decays = int(n_intervals*decay_per_interval)


def initial_energy_spectra(isotope,color):
    # digregard interval info for plotting
    JANIS_energy_matrix_keV = JANIS_properties(
        n_decays, isotope)[1]/1000  # convert to keV
    # flatten and remove zeros; interval dismissed
    JANIS_E_keV = np.ndarray.flatten(
        np.array(JANIS_energy_matrix_keV[JANIS_energy_matrix_keV != 0]))
    ### histogram ###
    # the bin size is based on smallest energy difference
    E_unique_ordered = np.sort(np.unique(JANIS_E_keV))
    E_diff = E_unique_ordered[1:] - E_unique_ordered[:-1]
    res = min(E_diff) # resolution or bin step size
    left_bound = min(JANIS_E_keV - res)
    right_bound = max(JANIS_E_keV + res)
    n_bins = int((right_bound - left_bound)/res) # range/res
    n, b = np.histogram(JANIS_E_keV, bins=n_bins, range=(left_bound, right_bound)) # generate spectra/bins
    # normalize counts
    norm_n = n/len(JANIS_E_keV)
    # bin midpoints
    b_mid = (b[:-1] + b[1:])/2
    # labels
    label1 = 'JANIS ' + isotope
    label2 = 'BRICC ' + isotope
    # final norm hist
    plt.hist(b_mid, range=(left_bound, right_bound),
             bins = n_bins, weights=norm_n, log=False, label=label1, color=color)
    plt.ylabel('Intensity')
    plt.xlabel('Auger Electron Energies (keV)')
    plt.legend()
    plt.tight_layout()
    #### need to add number of decays to graph
    

def spectra(isotope, LB,RB, color):
    # digregard interval info for plotting
    JANIS_energy_matrix_keV = JANIS_properties(
        n_decays, isotope)[1]/1000  # convert to keV
    # flatten and remove zeros; interval dismissed
    JANIS_E_keV = np.ndarray.flatten(
        np.array(JANIS_energy_matrix_keV[JANIS_energy_matrix_keV != 0]))
    ### histogram ###
    # the bin size is based on smallest energy difference
    E_unique_ordered = np.sort(np.unique(JANIS_E_keV))
    E_diff = E_unique_ordered[1:] - E_unique_ordered[:-1]
    res = min(E_diff)  # resolution or bin step size
    left_bound = min(JANIS_E_keV - res)
    right_bound = max(JANIS_E_keV + res)
    # print(left_bound,right_bound)
    n_bins = int((right_bound - left_bound)/res)  # range/res
    n,b = np.histogram(JANIS_E_keV,n_bins)
    b_mid = (b[:-1] + b[1:])/2
    norm_n = n/len(JANIS_E_keV)
    # # labels
    label1 = 'JANIS ' + isotope
    # label2 = 'BRICC ' + isotope
    plt.hist(b_mid, range=(left_bound, right_bound),
             bins=n_bins, weights=norm_n,log = True, label=label1,color=color)
    plt.xlim([LB,RB])
    # # final norm hist
    plt.ylabel('Intensity')
    plt.xlabel('Auger Electron Energies (keV)')
    plt.legend()
    plt.tight_layout()
    #### need to add number of decays to graph


#### select individually ####
# initial_energy_spectra('Pt195m','blue')
# initial_energy_spectra('HG197','red')
# initial_energy_spectra('Sb119','purple')
# initial_energy_spectra('Tb161','hotpink')

# spectra('I125',0,5,'blue')
# plt.show()
