
import matplotlib.pylab as plt
import numpy as np
from Energy_Data import auger_event_generator
from Energy_Data import get_auger_data
from scipy import signal

print('------------------------------------')
print('NEW RUN') 
print('------------------------------------')

'''
Parameters
'''

# options: 'beta', 'alpha',...
isotope = 'Pt195m'
# options: 'parallel plates', 'cylinder', 'box',...
geometry = 'PP'
# options: 'Si', 'Al'
material = 'Si'
# decay per interval
decay_per_interval = 10000
interval = 2
# number of particles
n_decays = int(interval*decay_per_interval)
# distance of souce from plate 1
source_displacement = 0.005e-3
# resolutions (eV)
res = [5,50,100]
nbins = interval*decay_per_interval*5
max_range = 2

'''
Materials
'''

def Si():
    k = 0.77
    J = 8.1517 # ionization potential (eV)
    Z = 14 # atomic number
    A = 18.085 # atomic weight amu
    p = 2.33 # density g/cc
    return [k,J,Z,A,p]

'''
Plate Geometry
'''

def parallel(source_displacement): 
    height = 40e-3 # width (x) 40mm
    width = 40e-3 # height (y)  40mm
    diameter = 0.4e-3 # 0.4 mm   # detector diameter (z)
    thickness = source_displacement*2 #thickness of the plate 0.02mm
    distance_source_to_inner_plate2 = diameter + source_displacement
    return [width,height,diameter,distance_source_to_inner_plate2,thickness]

'''
Isotope/Decay
'''
# Beta particle/fast electron

def Pt195m(isotope, n_decays):
    mass = 1    # mass is expressed in electron mass
    E_matrix = auger_event_generator(isotope,n_decays)
    n_electrons = np.shape(E_matrix)[1]
    return [mass,E_matrix,n_electrons]

'''
checking functions
'''
# fucntion that turns all non-zeros values to ones for checking
def tagged(v):
    tag = v != np.zeros(np.shape(v))
    v1 = v.copy()
    v1[tag] = 1
    return v1

def check(v0,v1,v2,v3,v4):
    t0 = tagged(v0)[0]
    t1 = tagged(v1)[0]
    t2 = tagged(v2)[0]
    t3 = tagged(v3)[0]
    t4 = tagged(v4)[0]
    c1 = t1 + t3 - t0
    c2 = t2 + t4 + t3 - t0
    c3 = t1 - t4 - t2
    sm = sum(c1 + c2 +c3)
    return sm == 0

'''
Decay Vectors
'''

def initialPositions(n_decays,n_electrons,source_displacement,geometry):
    if geometry == 'PP':
        width = parallel(source_displacement)[0] # width (x)
        height = parallel(source_displacement)[1] # height (y) 
        # note: source is always located at z = 0
    # randomize x,y,z,phi,theta matrices n = n_electrons, m = n_decays
    ones = np.ones([n_decays, n_electrons])
    x = (np.random.rand(n_decays,n_electrons))-(ones*0.5*width)
    y = (np.random.rand(n_decays, n_electrons))-(ones*0.5*height)
    z = np.zeros([n_decays,n_electrons])
    phi = (np.random.rand(n_decays,n_electrons))*np.pi         # test angle distribution
    theta = (np.random.rand(n_decays,n_electrons))*2*np.pi
    # combine into large position and momentum matrix
    source_vectors = np.vstack([[x,y,z]])
    momentum_direction = np.vstack([[theta,phi]])
    return [source_vectors,momentum_direction]

def projection(r0,dir,v_mag):
    # source vector components
    x0 = r0[0]
    y0 = r0[1]  
    z0 = r0[2]
    # momentum directions
    theta = dir[0]
    phi = dir[1]
    n_electrons = len(x0[0])
    # initiate vector component arrays
    x,y,z = [[0]*n_electrons]*n_decays,[[0]*n_electrons]*n_decays,[[0]*n_electrons]*n_decays
    # projection to length v_mag
    y += x0 + v_mag*np.cos(phi)*np.tan(theta)
    x += y0 + v_mag*np.sin(phi)*np.tan(theta)
    z += z0 + v_mag
    # concatneate
    r = np.vstack([[x],[y],[z]])
    return r

def finalPositions(geometry,source_vectors,momentum_direction): # positions on plate one outer side
    if geometry == 'PP':
        width = parallel(source_displacement)[0] # width (x)
        height = parallel(source_displacement)[1] # height (y)
        distance_source_to_inner_plate2 = parallel(source_displacement)[3] # detector diameter (z)
        thickness = parallel(source_displacement)[4] # plate thickness
        phi = momentum_direction[1]
    # vector lengths
    distance_source_to_outer_plate1 = -(thickness-source_displacement) #points down
    distance_source_to_inner_plate1 = thickness-source_displacement # points up
    # final vectors
    inner_plate1 = projection(source_vectors,momentum_direction,distance_source_to_inner_plate1)
    outer_plate1 = projection(source_vectors,momentum_direction,distance_source_to_outer_plate1)
    inner_plate2 = projection(source_vectors ,momentum_direction,distance_source_to_inner_plate2)
    missed_plate2 = projection(source_vectors ,momentum_direction,distance_source_to_inner_plate2)
    # sorting out conditions
    up  = np.vstack([[(0 < phi) & (phi < np.pi/2)]*3])
    down = ~up
    inside  =  np.vstack([[(abs(inner_plate2[0]) < height/2) &  (abs(inner_plate2[1]) < width/2) ]*3])
    outside = ~inside
    # apply conditions to vectors
    inner_plate1[down] = 0  # remove decays going down
    outer_plate1[up] = 0  # remove decays going up
    inner_plate2[down] = 0 # remove decays going down
    inner_plate2[outside] = 0 # remove decays not on plate2
    missed_plate2[down] = 0 # remove decays going down
    missed_plate2[inside] = 0 # remove decays on plate2
    #check
    ch = check(source_vectors[0],inner_plate1[0],inner_plate2[0],outer_plate1[0],missed_plate2[0])
    if ch == True:
        print('particles were conserved')
    else:
        print('Error: particles were not conserved')
    return [[inner_plate1],[inner_plate2],[outer_plate1],[missed_plate2]]

'''
Energies
'''
def displacement(r0,r1,E):
    n_electrons = np.shape(E)[1]
    rf = np.vstack([[r1[0]-r0[0]],[r1[1]-r0[1]],[r1[2]-r0[2]]])
    rf_matrix = np.array(np.sqrt(rf[0]**2+rf[1]**2+rf[2]**2))
    E = E.copy() 
    tag = E == np.zeros(np.shape(E))
    rf_matrix[tag] = 0
    r_matrix = np.reshape(rf_matrix,(n_decays,n_electrons))
    return r_matrix


def SpecificEnergyLoss(E,material_values):
    # material Values
    k = material_values[0] # material dep value
    J = material_values[1] # mean ionization potential
    Z = material_values[2] # atomic number
    A = material_values[3] # atomic weight
    p = material_values[4] # density (g/cc)
    # copy E and make 0 = -1 to relove /0 error
    tag = E == 0
    E = np.copy(E)
    E[tag] = -1
    # compute the first part of the stopping power
    v1 = 785*((p*Z)/(A*E))
    v1[tag] = 0
    SP = v1*np.log((1.166*(E+k*J))/J)
    return SP

def energyLost(E,SP,R):
    v1 = np.asarray(SP).reshape(-1) * np.asarray(R).reshape(-1)
    E_lost = np.reshape(v1,[np.shape(E)[0],np.shape(E)[1]])         # SP*R
    E_final = E - E_lost     #E_initial - SP*R
    E = np.matrix(E.copy()) 
    tag = E == np.zeros(np.shape(E))
    E_final[tag] = 0
    return E_final

def energyDetected(n_decays, E_matrix, material_values,IP,FP,mass):
    # decay vectors FP = [inner_plate1,inner_plate2,outer_plate1,missed_plate2]
    source_vectors = IP[0]
    inner_plate1 = FP[0].copy()
    inner_plate2  = FP[1].copy()
    outer_plate1 = FP[2].copy()
    missed_plate2 = FP[3].copy()
    # detected conditions using positions
    n_electrons = len(FP[0][0][0])
    in_P1_E = inner_plate1[0] == np.zeros([n_decays,n_electrons])
    in_P2_E = inner_plate2[0] == np.zeros([n_decays,n_electrons])
    out_P1_E = outer_plate1[0] == np.zeros([n_decays,n_electrons])
    miss_P2_E = missed_plate2[0] == np.zeros([n_decays,n_electrons])
    # create new arrays for Energies
    E0 = E_matrix.copy()
    E1 = E_matrix.copy()
    E2 = E_matrix.copy()
    E3 = E_matrix.copy()
    E4 = E_matrix.copy()
    # sort out energy based on conditions
    E1[in_P1_E] = 0
    E2[in_P2_E] = 0
    E3[out_P1_E] = 0
    E4[miss_P2_E] = 0
    # check energies
    ch = check(E_matrix[0],E1[0],E2[0],E3[0],E4[0])
    x = ch == [[True]*n_electrons]
    if x.all() == True:
        print('Energies were conserved')
    else:
        print('Error: Energies were not conserved')
    return [[E0],[E1],[E2],[E3],[E4]]

# redo take range from 0 to 10 for smaller range

def HistE(Sm,fig,s,n):
    tabulated_energies_eV, tabulated_intensities = get_auger_data(isotope)
    tabulated_energies_eV= np.array(tabulated_energies_eV)*(1/1000)
    E = Sm[Sm != 0]/1000
    fig.set_size_inches(9,6)
    limits = np.arange(0,n_decays+1,decay_per_interval)
    plt.plot(tabulated_energies_eV,tabulated_intensities,'x',markersize = 10,color = 'orange',label = 'Expected Spectra')
    E = np.ndarray.flatten(E)
    n_spec, nbin = np.histogram(E, bins = 1000,range = (0,240))
    exc = n_spec != 0
    plt.errorbar(nbin[1:][exc],(n_spec/max(n_spec)*331)[exc],yerr = 331*np.sqrt(max(n_spec)**2 * n_spec + n_spec**2 * max(n_spec))[exc]/(max(n_spec**2)), fmt = 'o',color = 'blue',label = 'Simulated Spectra')
    plt.yscale('log')
    plt.ylabel('Intensity')
    plt.xlabel('Final Particle Energy (keV)')
    plt.title(f'Energies Detected on Plate {n}')
    plt.legend(loc = 'upper right')
    plt.figtext(0.05,0.0005,s)

def Hist2D(E1,E2):
    #print(E1[:10], E2[:10])
    # ratio = len(E2)/len(E1)
    # nbins2 = int(240*ratio)
    # n1, b1 = np.histogram(E1, bins = 100,range = (0,240))
    # n2, b2 = np.histogram(E2, bins = 100,range = (0,240))
    # bins1 = (b1[1:] + b1[:-1])/2
    # bins2 = (b2[1:] + b2[:-1])/2
    plt.hist(E1, bins = 120,range = (0,20), color = 'blue', label = 'Plate 1')
    plt.hist(E2, bins = 120, range = (0,20), color='red', label = 'Plate 2')
    plt.xlabel('Deposited Energy keV')
    plt.ylabel('Intensity')
    plt.legend()
    plt.figure()
    plt.hist2d(E1, E2, bins = (100,100), range = ((0,20),(0,20)), cmap='Greens')
    plt.xlabel("Plate 1 (- what doesn't make it to Plate2)")
    plt.ylabel('Plate 2')
    #plt.scatter(E1, E2)

def stoppingPowerGraph(E0,SP_exp,dimensions,n_electrons,material_values):
    # build theoretical data
    E0 = np.asarray(np.matrix.flatten(np.matrix(E0))).reshape(-1)
    SP_exp = np.asarray(np.matrix.flatten(np.matrix(SP_exp))).reshape(-1)
    E0_theory = np.arange(0,240000,1)
    #E0_theory = np.linspace(min(E0),max(E0),(int(n_electrons*n_decays*10)))
    r_max = np.sqrt(dimensions[0]**2+dimensions[1]**2+dimensions[2]**2)
    print('r = ', r_max)
    SP_theory = SpecificEnergyLoss(E0_theory,material_values)
    E_lost_theory = E0_theory - SP_theory*r_max
    # caluclate experimental
    E_lost_exp = E0 - SP_exp*r_max
    # get up subplots
    fig = plt.figure(0)
    fig.set_size_inches(12,6)
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    # remove zeros
    E = np.copy(E0)
    tag = E != 0
    E0 = E[tag]
    E_lost_exp = E_lost_exp[tag]
    SP_exp = SP_exp[tag]
    # theory graph
    ax1.plot(E0_theory/1e3,SP_theory,'--',color ='black',label = 'Model')
    ax2.plot(E0_theory/1e3,E_lost_theory/1e3,'--',color ='black',label = 'Model')
    # experimental graph
    ax1.plot(E0/1e3,SP_exp,'o',color = 'red',label = 'Simulated')
    ax2.plot(E0/1e3,E_lost_exp/1e3,'o',color = 'red',label = 'Simulated')
    # titles and labels
    ax1.set_title('Specific Energy Loss vs Initial Particle Energy')
    ax1.set_ylabel('Specific Energy Loss (MeV/mm)')
    ax1.set_xlabel('Initial Particle Energy (keV)')
    ax2.set_title('Final vs Initial Electron Energy')
    ax2.set_xlabel('Electron Initial Energy (keV)')
    ax2.set_ylabel('Final Electron Energy (keV)')
    ax1.legend(loc = 'upper right')
    ax2.legend(loc = 'upper left')

def guassian(x,a,u,s):
    return a*abs(x[1] - x[0])*np.exp(-((x-u)**2)/(2*s**2))/(np.sqrt(2*np.pi*s**2))

def convolve(E,res,nbins,max_range,fig1,s):
    E = np.ndarray.flatten(np.array(E[E!=0]))/1000
    # Experimental Error
    noise1 = np.random.normal(0,res[0],len(E))/1000
    noise2 = np.random.normal(0,res[1],len(E))/1000
    noise3 = np.random.normal(0,res[2],len(E))/1000
    Em1 = E + noise1
    Em2 = E + noise2
    Em3 = E + noise3
    # data for normalization
    y, x = np.histogram(E, bins = nbins,range = (0,max_range))
    ym1, xm1 = np.histogram(Em1, bins = nbins,range = (0,max_range))
    ym2, xm2 = np.histogram(Em2, bins = nbins,range = (0,max_range))
    ym3, xm3 = np.histogram(Em3, bins = nbins,range = (0,max_range))
    x =  np.asarray((x[1:] + x[:-1])/2).reshape(-1)
    y =  np.asarray(y).reshape(-1)
    xm1 =  np.asarray((xm1[1:] + xm1[:-1])/2).reshape(-1)
    ym1=  np.asarray(ym1).reshape(-1)
    xm2 =  np.asarray((xm2[1:] + xm2[:-1])/2).reshape(-1)
    ym2=  np.asarray(ym2).reshape(-1)
    xm3 =  np.asarray((xm3[1:] + xm3[:-1])/2).reshape(-1)
    ym3=  np.asarray(ym3).reshape(-1)
    # create filters at 3 given resolutions
    filt1 = guassian(np.linspace(-res[0]*10,res[0]*10+1),1,0,res[0]) # should be in eV because of change in binning, units of gausian match histogram units 
    filt2 = guassian(np.linspace(-res[1]*10,res[1]*10+1),1,0,res[1])
    filt3 = guassian(np.linspace(-res[2]*10,res[2]*10+1),1,0,res[2])
    # create figure
    fig1.set_size_inches(12,4.5)
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133)
    # figure 1 graphs with noise
    norm1 = 1
    norm2 = 1
    norm3 = 1
    ax1.plot(xm1,signal.fftconvolve(ym1,filt1,'same')*norm1,'-',color = 'red')
    ax2.plot(xm2,signal.fftconvolve(ym2,filt2,'same')*norm2,'-',color = 'red')
    ax3.plot(xm3,signal.fftconvolve(ym3,filt3,'same')*norm3,'-',color = 'red',label = 'detector response with noise')
    # titles and axes
    ax1.set_title(f'Resolution = {res[0]} eV')
    ax1.set_xlabel('Final Particle Energy (keV)')
    ax1.set_ylabel('Intensity')
    ax1.set_xlim([0,2])
    ax2.set_title(f'Resolution = {res[1]} eV')
    ax2.set_xlabel('Final Particle Energy (keV)')
    ax2.set_xlim([0,2])
    ax3.set_title(f'Resolution = {res[2]} eV')
    ax3.set_xlabel('Final Particle Energy (keV)')
    ax3.set_xlim([0,2])
    plt.figtext(0.05,0.0005,s)

# if same binning for all graphs put in the units of y axis enter in 1/binwidth or normailzation (divide by binwidth)
### check a plot of the diffference between simulated with nois and simulated, should be random if correct

'''
Results
'''
def detector(n_decays,source_displacement,geometry,isotope,material,res,nbins,max_range):
    if isotope == 'Pt195m':
        decays = Pt195m(isotope, n_decays)
        n_electrons = decays[2]
        E_matrix = decays[1]
        mass = decays[0]
    if geometry == 'PP':
        PP = parallel(source_displacement)
        dimensions = [PP[0],PP[1],PP[2],PP[3],PP[4]]
        IP = initialPositions(n_decays,n_electrons,source_displacement,geometry)
        source_vectors = IP[0]
        momentum_directions = IP[1]
        FP = np.vstack(finalPositions(geometry,source_vectors,momentum_directions))
    if material == 'Si':
        material_values = Si()
    E_array = np.vstack(energyDetected(n_decays, E_matrix, material_values,IP,FP,mass))
    # Specific Energy Loss
    E0 = E_array[0]
    E1 = E_array[1]
    E2 = E_array[2]
    E3 = E_array[3]
    E4 = E_array[4]
    # E1 Energy loss
    SP1 = SpecificEnergyLoss(E1,material_values)
    R1 = displacement(source_vectors,FP[0],E1)
    E1_lost = energyLost(E1,SP1,R1)
    # E2 Energy Loss
    tag = E2 == 0
    E2_lost = np.copy(E1_lost)
    E2_lost[tag] = 0 
    # E3 Energy Loss
    SP3 = SpecificEnergyLoss(E3,material_values)
    R3 = displacement(source_vectors,FP[2],E1)
    E3_lost = energyLost(E3,SP3,R3)
    # E4 Energy Loss
    tag = E4 == 0
    E4_lost = np.copy(E1_lost)
    E4_lost[tag] = 0 
    # total values
    E_lost_tot_0 = np.array(E1_lost + E3_lost)
    E_lost_tot = E_lost_tot_0[E_lost_tot_0 != 0]
    E_array = [E0,E1_lost[E1_lost!=0],E2_lost[E2_lost!=0],E3_lost[E3_lost!=0],E4_lost[E4_lost!=0]] 
    # SP_total = np.array(SP1 +SP3)
    # stoppingPowerGraph(E0,SP_total,dimensions,n_electrons,material_values)
    print(len(E2_lost[E2!=0])/(len(E0[E0!=0]))*100, 'percent of all the electrons were detected on plate 2.')
    print('number E2 = ', len(E2_lost[E2!=0]))
    print('total electrons = ', len(E0[E0!=0]))
    s = f'number of auger electrons detected on plate 2 = {round(len(E2_lost[E2!=0]),1)}, {round(len(E2_lost[E2!=0])/(len(E0[E0!=0]))*100,1)} percent of all electrons , material = ' + material + ', isotope = ' + isotope
    s1 = f'number of auger electrons detected on plate 1 = {round(len(E_lost_tot),1)}, {round(len(E_lost_tot)/(len(E0[E0!=0]))*100,1)} percent of all electrons , material = ' + material + ', isotope = ' + isotope
    detector2 = np.asarray(E_array[2]).reshape(-1)
    detector1 = np.asarray(E_lost_tot).reshape(-1)
    # 2D Histogram
    E0_flat = np.asarray(E0).reshape(-1)
    E2_flat = np.asarray(E2).reshape(-1)
    E1_flat = np.asarray(E1).reshape(-1)
    Etot_flat = np.asarray(E_lost_tot_0).reshape(-1)
    E2_lost_flat = np.asarray(E2_lost).reshape(-1)
    mask = E1_flat != 0
    plate1 = E0_flat[mask]
    plate2 = E2_flat[mask]
    #Spectra Graphs
    fig2 = plt.figure(2)
    print('d1 = ',detector1[:10],'p1 =',plate1[:10])
    HistE(plate1, fig2, s,1)
    fig3 = plt.figure(3)
    HistE(plate2, fig3, s1,2)
    fig6 = plt.figure(6)
    Hist2D(plate1/1000,plate2/1000)
    cb = plt.colorbar()
    cb.set_label('Intensity')
    # Detector resonse Graphs
    fig4 = plt.figure(4)
    convolve(plate2,res,nbins,max_range,fig4,s)
    fig5 = plt.figure(5)
    convolve(plate1,res,nbins,max_range,fig5,s1)
    plt.tight_layout()
    plt.show()
    return [IP,FP,E_array]

Results = detector(n_decays,source_displacement,geometry,isotope,material,res,nbins,max_range)

print('end of code')











