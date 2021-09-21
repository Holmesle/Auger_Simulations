import pandas as pd
import numpy as np



def get_auger_data(isotope):
    tabulated_energies_eV = []
    tabulated_intensities = []
    data_filename = "/Users/holm423/OneDrive - PNNL/Desktop/SULI Intern Project/Code/Janis Energy Data/janis_{}.csv".format(isotope)
    with open(data_filename,'r') as myfile:
        d = dict(pd.read_csv(myfile, delimiter=';'))
        for i in range(len(d['E '])):
            tabulated_energies_eV.append(float(str(d['E '][i]).split(" eV")[0]))
            tabulated_intensities.append(float(str(d[' I '][i]).split(" %")[0]))  
    return tabulated_energies_eV, tabulated_intensities


def auger_event_generator(isotope,n_decays):
    tabulated_energies_eV, tabulated_intensities = get_auger_data(isotope)
    intensity_matrix = np.matrix(np.random.poisson(lam = (np.array(tabulated_intensities[:])/100),size = [n_decays,len(tabulated_intensities)]))
    decay_energies_eV = []
    sum_array = intensity_matrix.sum(axis=1)
    n_electrons = np.copy(np.matrix.max(sum_array))
    energy_matrix = []
    n_intensities = np.shape(intensity_matrix)[1]
    for i in range(0,n_decays):
        for j in range(0,n_intensities):
            if intensity_matrix[i,j] != 0:
                decay_energies_eV += [tabulated_energies_eV[j]]*intensity_matrix[i,j]
        L = n_electrons - len(decay_energies_eV)
        if L != 0:
            decay_energies_eV += [0]*L
        energy_matrix += [decay_energies_eV]
        decay_energies_eV =[]
    energy_matrix = np.matrix(energy_matrix)
    return energy_matrix
