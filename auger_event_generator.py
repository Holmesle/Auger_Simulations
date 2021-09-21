import numpy as np
import csv
#import ROOT as root



def get_auger_data(isotope='Ga67'):

    # create lists to hold the energies and intensities read from the data files.
    tabulated_energies_eV = []
    tabulated_intensities = []
    
    # a list to contain the energies of all the electrons emitted in a single decay.
    # This is the return value of this function.
    decay_energies_eV = []

    # Get the contents of each row of the CSV as an OrderedDict with keys given by the first row.
    # Parse the rows to extract the energy and intensity values into lists
    data_filename = '/Users/holm423/OneDrive - PNNL/Desktop/Research/janis_{}.csv'.format(isotope)
    with open(data_filename, newline='') as csvfile:

        dict_reader = csv.DictReader(csvfile, delimiter=";")

        for row in dict_reader:

            # note the importance of the leading and trailing spaces below. JANIS wites
            # csv as 'xxxx ; yyyy ; zzzz'    
            tabulated_energies_eV.append(float(row["E "].split(" ")[0]))
            tabulated_intensities.append(0.01*float(row[" I "].split(" ")[1]))

    return tabulated_energies_eV, tabulated_intensities


def auger_event_generator(tabulated_energies_eV, tabulated_intensities):
    
    # loop over every electron in the list, and decide how many of them we get in this decay.
    # Assume a model of Auger/CK/CE decay where each electron is emitted independently 
    # of all the others. That is probably a terrible model that should be refined. The number of each
    # electron emitted is given by drawing a random number from a poisson distribution with the 
    # tabulated intensity as the mean.
    
    decay_energies_eV = []
    
    # for every possible electron energy
    for i in range(len(tabulated_energies_eV)):

        # decide how many electrons of this energy
        n_electrons = np.random.poisson(tabulated_intensities[i])
        
        # Put those electron energies into the returned list
        for j in range(n_electrons):
            decay_energies_eV.append(tabulated_energies_eV[i])
            
            
    # return a list containing the energies of every electron emitted in a single decay.
    return decay_energies_eV



# Now test that model. 

# Choose an isotope
isotope = "Hg197"

# Get the data from file
tabulated_energies_eV, tabulated_intensities = get_auger_data(isotope)

# Create a ROOT histogram to contain the spectrum
#energy_spectrum = root.TH1F("energy_spectrum_{}".format(isotope), "{} Electron Energy Spectrum".format(isotope), 1000, 0., 6.)

# Loop over individual decay events
n_decays = 1000

for i in range(n_decays):
    
    # a list to contain the energies of all the electrons emitted in a single decay.
    decay_energies_eV = auger_event_generator(tabulated_energies_eV, tabulated_intensities)
    
    # loop over the electron energies returned by auger_event_generator and enter them into a histogram
    for j in range(len(decay_energies_eV)):
        energy_spectrum.Fill(np.log10(decay_energies_eV[j]))
        
    if (i%100 == 0):
        print(i, decay_energies_eV)
            
            
# draw it
#%jsroot on
#tcanvas = root.TCanvas("tcanvas", "Auger Energy Spectrum", 800, 600)
#energy_spectrum.Draw()
#tcanvas.Draw()
