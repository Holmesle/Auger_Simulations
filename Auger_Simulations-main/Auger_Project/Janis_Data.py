import pandas as pd
import numpy as np

'''
This module requires no module inputs but outputs an initial_energy_matrix from the the Janis csv
library to the isotope_source module.

     BRICC -> Isotope_Source -> Detector
Janis_Data -> 

'''

def get_auger_data(isotope):
    data_filename = r"C:\Users\holme\OneDrive\Desktop\python_code\Auger_Simulations-main\Janis\janis_{}.csv".format(
        isotope)
    with open(data_filename, 'r') as myfile:
        d = dict(pd.read_csv(myfile, delimiter=';'))
        E = np.array((d['E ']).str.split(' eV'))
        I = np.array((d[' I ']).str.split(' %'))
        tabulated_energies_eV = np.hstack(E)[0::2].astype(float)
        tabulated_intensities = np.hstack(I)[0::2].astype(float)
    return tabulated_energies_eV,tabulated_intensities

def initial_energies(isotope,n_decays):
    tabulated_energies_eV, tabulated_intensities = get_auger_data(isotope)
    # determine how many electrons of each energy are generated (n_decays,n_energies)
    intensity_matrix = np.array(np.random.poisson(lam = (np.array(tabulated_intensities[:])/100),size = [n_decays,len(tabulated_intensities)]))
    # find max number of electrons generated in each decay 
    # to determine final matrix size (n_decays,max n_electrons)
    sum_array = intensity_matrix.sum(axis=1)
    n_electrons = np.amax(sum_array)
    # 
    energy_matrix = []
    decay_energies_eV = []
    n_intensities = np.shape(intensity_matrix)[1]
    for i in range(0,n_decays):
        for j in range(0,n_intensities):
            if intensity_matrix[i,j] != 0:
                # each energy value will appear as many times as the intiensity value
                decay_energies_eV += [tabulated_energies_eV[j]]*intensity_matrix[i,j]
        L = n_electrons - len(decay_energies_eV)
        if L != 0:
            decay_energies_eV += [0]*L
        energy_matrix += [decay_energies_eV]
        decay_energies_eV =[]
    energy_matrix = np.array(energy_matrix)
    energy_matrix = np.hstack(energy_matrix)
    energy_matrix = np.reshape(energy_matrix,[n_decays,n_electrons])
    return energy_matrix


# print(initial_energies('Pt195m', 10))

