import numpy as np
import matplotlib.pylab as plt
from Energy_Data import get_auger_data

isotopes = ['Ga67','Hg197m','I125','In111','Pt195m']

print('start ----------------------------')


def Spectra(isotope,max_range,nbins,color):
    tabulated_energies_eV, tabulated_intensities = get_auger_data(isotope)
    E = []
    for i in range(len(tabulated_energies_eV)):
        if round(tabulated_intensities[i]) != 0:
            #print(round(tabulated_intensities[i]*1000))
            E += [tabulated_energies_eV[i]/1000]*round(tabulated_intensities[i]*100)
    E = np.array(E)
    E = np.ndarray.flatten(np.array(E))
    n, b = np.histogram(E, bins = nbins,range = (0,max_range))
    plt.hist((b[:-1] + b[1:])/2, range = (0,max_range), weights = n ,bins = nbins,log = True, label = isotope,color = color)
    #plt.plot(np.array(tabulated_energies_eV)/1000, tabulated_intensities,'x')


fig = plt.figure()
Spectra('Pt195m',240,480,'blue')
Spectra('Hg197',900,1800,'red')
Spectra('I123',600,1200,'orange')
Spectra('Sb119',280,280*2,'purple')
Spectra('Tb161',550,550*2,'hotpink')
Spectra('Tc99m',330,330*2,'gold')
plt.ylabel('Intensity')
plt.xlabel('Auger electron Energies (keV)')
plt.xlim([0,100])
#plt.ylim([0,5000])
plt.title('Tabulated Isotope Intensities (below 100 keV)')
plt.legend(loc = 'upper right')
print('end ----------------------------')
plt.show()
