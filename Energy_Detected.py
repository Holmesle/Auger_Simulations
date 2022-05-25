import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Geometry as GM
import Isotope_Source as IS
import Materials as MT
import Particle_Path as PTH
import Specific_Energy_Loss as SEL
import Parallel_Plate as PP 
import csv

# Energy and Particles
n_decays = 1
isotope = 'I125'
m, E = IS. JANIS_properties(n_decays, isotope)
n_rows = len(E[0,:])
n_particles = len(E.reshape(n_rows*n_decays)[E.reshape(n_rows*n_decays)!=0])
# Materials j
material = 'Si'
material_values = MT.material_dict(material)
# Plate Dimensions
height = 40e-3  # width (x) 40mm
width = 40e-3  # height (y)  40mm
source_to_inner_surface = 0.005e-3
source_to_outer_surface = 1e-1
thickness = source_to_inner_surface  + source_to_outer_surface
diameter = 0.4e-3  # 0.4 mm   # detector diameter (z)

# Geometry
dim, pos, angle, pos_matrices, angle_matrices = GM.Parallel_Plate(E, height, width, source_to_inner_surface, source_to_outer_surface, diameter)

'''
Useful and Frequent Functions 
'''
# function that removes zeros and make array 1 dimensional
def remove_zeros(v,shape):
    if v.ndim == 3:
        x0 = v[0].reshape(shape)
        y0 = v[1].reshape(shape)
        z0 = v[2].reshape(shape)
        sum0 = x0+y0+z0
        x0 = x0[sum0 != 0]
        y0 = y0[sum0 != 0]
        z0 = z0[sum0 != 0]
        return np.vstack([x0,y0,z0])
    else:
        v = v.reshape(shape)
        return v[v!=0]

# function that turns all non-zero entries into 1 
# and reshapes vector into the given shape
def make_filter(v, shape):
    v= v.copy()
    if v.ndim > 2:
        v = sum(v)
    v[v==0]=0
    v[v!=0]=1
    v = v.reshape(shape)
    return v

def write(Type_E,Es):
    file = r"C:\Users\holme\OneDrive\Desktop\python_code\Auger_Simulations-main\Auger_Project\check_Energy.csv"
    f = open(file,'w')
    writer = csv.writer(f)
    for i in range(len(Type_E)):
        row = [Type_E[i]] + [str(j) for j in Es[i]]
        writer.writerow(row)
    f.close()

proj, E_array = PTH.Paths(E,n_decays,n_rows,pos_matrices,angle_matrices,dim)

proj = proj.reshape(5,3,n_decays,n_rows)
proj1 = proj[0]
proj2 = proj[1]
proj3 = proj[2]
proj4 = proj[3]
proj5 = proj[4]
proj1_reduced = proj1[proj1!=0]
proj2_reduced = proj1[proj2!=0]
proj3_reduced = proj1[proj3!=0]
proj4_reduced = proj1[proj4!=0]
proj5_reduced = proj1[proj5!=0]

proj_plot = [proj1_reduced,proj2_reduced,proj3_reduced,proj4_reduced,proj5_reduced]

def Initial_Positions(proj):
    # unpack proj
    proj = proj.reshape(5,3,n_decays,n_rows)
    proj1 = proj[0]
    proj2 = proj[1]
    proj3 = proj[2]
    proj4 = proj[3]
    proj5 = proj[4]
    # find initial positions of each population
    pos1 = pos_matrices.copy().reshape(3,n_decays,n_rows)
    pos2 = pos_matrices.copy().reshape(3,n_decays,n_rows)
    pos3 = pos_matrices.copy().reshape(3,n_decays,n_rows)
    pos4 = pos_matrices.copy().reshape(3,n_decays,n_rows)
    pos5 = pos_matrices.copy().reshape(3,n_decays,n_rows)
    pos1[proj1==0] = 0
    pos2[proj2==0] = 0
    pos3[proj3==0] = 0
    pos4[proj4==0] = 0
    pos5[proj5==0] = 0
    initial_pos = [pos1,pos2,pos3,pos4,pos5]
    return initial_pos

def Final_Positions(proj,E_array):
     # unpack proj
    proj = proj.reshape(5,3,n_decays,n_rows)
    proj1 = proj[0]
    proj2 = proj[1]
    proj3 = proj[2]
    proj4 = proj[3]
    proj5 = proj[4]
    # unpack E
    E1 = E_array[0:n_decays]
    E2 = E_array[n_decays:n_decays*2]
    E3 = E_array[n_decays*2:n_decays*3]
    E4 = E_array[n_decays*3:n_decays*4]
    E5 = E_array[n_decays*4:n_decays*5]
    # max distance calculated from SEL
    s1 = SEL.max_s(E1,n_rows,n_decays)*1e-7 #convert to milimeters from angstoms
    s2 = SEL.max_s(E2,n_rows,n_decays)*1e-7 
    s3 = SEL.max_s(E3,n_rows,n_decays)*1e-7 
    s4 = SEL.max_s(E4,n_rows,n_decays)*1e-7 
    s5 = SEL.max_s(E5,n_rows,n_decays)*1e-7 
    # find vectors using s matrix
    s1_v = PTH.projection(pos_matrices,angle_matrices,-s1,n_decays).reshape(3,n_decays,n_rows)
    s2_v = PTH.projection(pos_matrices,angle_matrices,s2,n_decays).reshape(3,n_decays,n_rows)
    s3_v = PTH.projection(pos_matrices,angle_matrices,s3,n_decays).reshape(3,n_decays,n_rows)
    s4_v = PTH.projection(pos_matrices,angle_matrices,s4,n_decays).reshape(3,n_decays,n_rows)
    s5_v = PTH.projection(pos_matrices,angle_matrices,s5,n_decays).reshape(3,n_decays,n_rows)
    # get rid of values not in pop
    s1_v[proj1==0]=0
    s2_v[proj2==0]=0
    s3_v[proj3==0]=0
    s4_v[proj4==0]=0
    s5_v[proj5==0]=0
    # filters to find particles in plate 
    f1_inx = np.array([(s1_v[0] <= dim[1]/2) & (s1_v[0] >= -dim[1]/2)]*3)
    f1_iny = np.array([(s1_v[1] <= dim[0]/2) & (s1_v[1] >= -dim[0]/2)]*3)
    f1_ext = np.array([s1_v[2] >= -dim[3]]*3)
    f2_inx = np.array([(s2_v[0] <= dim[1]/2) & (s2_v[0] >= -dim[1]/2)]*3)
    f2_iny = np.array([(s2_v[0] <= dim[1]/2) & (s2_v[0] >= -dim[1]/2)]*3)
    f2_int = np.array([s2_v[2] <= dim[2]]*3)
    f4_inx = np.array([(s4_v[0] <= dim[1]/2) & (s4_v[0] >= -dim[1]/2)]*3)
    f4_iny = np.array([(s4_v[0] <= dim[1]/2) & (s4_v[0] >= -dim[1]/2)]*3)
    f4_int = np.array([s4_v[2] <= dim[2]+dim[5]+dim[4]]*3)
    # final positions in plate 1
    fp1 = proj1.copy()
    fp2 = proj2.copy()
    fp3 = proj3.copy()
    fp4 = proj4.copy()
    fp5 = proj5.copy()
    print('s1_v <',dim[1]/2)
    print(s1_v[0])
    # everything in side plate make zero in s matrices
    s1_v[~f1_inx]=0 
    s1_v[~f1_iny]=0 
    s1_v[~f1_ext]=0 
    print(f1_inx)
    print('s1_v <',dim[1]/2)
    print(s1_v[0])
    s2_v[~f2_inx]=0 
    s2_v[~f2_iny]=0 
    s2_v[~f2_int]=0 
    s4_v[~f4_inx]=0 
    s4_v[~f4_iny]=0 
    s4_v[~f4_int]=0 
    s5_v[~f4_inx]=0 
    s5_v[~f4_iny]=0 
    s5_v[f4_int]=0 
    # proj1 replace non-zero entries with sn_v
    fp1[s1_v!=0] = s1_v[s1_v!=0]
    fp2[s2_v!=0] = s2_v[s2_v!=0]
    fp4[s4_v!=0] = s4_v[s4_v!=0]
    fp5[s5_v!=0] = s5_v[s5_v!=0]
    # high enough energy to pass plates surface add extra distance
    fp1z = fp1[2]
    fp2z = fp2[2]
    fp5z = fp5[2]
    fp3z = fp3[2]
    fp1z[~f1_ext[2]] = fp1z[~f1_ext[2]]  - dim[2]
    fp2z[~f2_int[2]] = fp2z[~f2_int[2]]  + dim[2]
    fp5z[fp5z!=0] = fp5z[fp5z!=0]  + dim[2]
    fp3z[fp3z!=0] = fp3z[fp3z!=0]  + dim[2]
    # remove pop3-5 form pop2
    # fp2[fp5!=0]=0
    # fp2[fp4!=0]=0
    # fp2[fp3!=0]=0
    # fp4 remaining in plate 2
    # fp4[fp5!=0]=0
    # some reason after line 189 a random zero appears no clue why
    fps = [fp1,fp2,fp3,fp4,fp5]
    return fps
    
def Conservation_of_Mass(E_array,fps):
    # unpack E
    E1, E2, E3, E4, E5 = E_array[0:n_decays], E_array[n_decays:n_decays*2], E_array[n_decays*2:n_decays*3], E_array[n_decays*3:n_decays*4], E_array[n_decays*4:n_decays*5]
    # unpack fps
    fp1, fp2, fp3, fp4, fp5 = fps
    # copy initial energies to remove placce holders
    E1_pop = make_filter(E1,n_decays*n_rows)
    E2_pop = make_filter(E2,n_decays*n_rows)
    E3_pop = make_filter(E3,n_decays*n_rows)
    E4_pop = make_filter(E4,n_decays*n_rows)
    E5_pop = make_filter(E5,n_decays*n_rows)
    # find partilces in each pop based on final position
    fp1_particles = make_filter(fp1,n_decays*n_rows)
    fp2_particles = make_filter(fp2,n_decays*n_rows)
    fp3_particles = make_filter(fp3,n_decays*n_rows)
    fp4_particles = make_filter(fp4,n_decays*n_rows)
    fp5_particles = make_filter(fp5,n_decays*n_rows)
    # check conservation of mass
    P1_particles = fp1_particles + fp2_particles
    P1_Energies = E1_pop + E2_pop
    P1_particles = P1_particles[P1_Energies!=0]
    n1 = remove_zeros(fp1_particles,(n_decays,n_rows)).size
    n2 = remove_zeros(fp2_particles,(n_decays,n_rows)).size
    n3 = remove_zeros(fp3_particles,(n_decays,n_rows)).size
    n4 = remove_zeros(fp4_particles,(n_decays,n_rows)).size
    n5 = remove_zeros(fp5_particles,(n_decays,n_rows)).size
    n_P1 = n1 + n2
    n_post_P1 = n3 + n4 + n5
    if n_P1 - n_particles != 0:
        print('Error: mass not conserved in P1')
        print('n total:', n_particles)
        print('n in P1:', n_P1)
    if n_post_P1 - n2 != 0:
        print('Error: mass not conserved in P2')
        print('n total:', n2)
        print('n after P1:', n_post_P1)
    n_P1 = n1 + n2
    n_P2 = n4 + n5
    n_escape = n3
    return n_P1,n_P2,n_escape


def Final_Energy(E_array,proj,fps, initial_pos):
    E1, E2, E3, E4, E5 = E_array[0:n_decays], E_array[n_decays:n_decays*2], E_array[n_decays*2:n_decays*3], E_array[n_decays*3:n_decays*4], E_array[n_decays*4:n_decays*5]
    fp1, fp2, fp3, fp4, fp5 = fps
    pos1, pos2, pos3, pos4, pos5 = initial_pos
    proj = proj.reshape(5,3,n_decays,n_rows)
    proj1 = proj[0]
    proj2 = proj[1]
    proj3 = proj[2]
    proj4 = proj[3]
    proj5 = proj[4]
    # find partilces in each pop based on final position
    fp1_particles = make_filter(fp1,n_decays*n_rows)
    fp2_particles = make_filter(fp2,n_decays*n_rows)
    fp3_particles = make_filter(fp3,n_decays*n_rows)
    fp4_particles = make_filter(fp4,n_decays*n_rows)
    fp5_particles = make_filter(fp5,n_decays*n_rows)
    # sort particles based on final position
    E1_filter = make_filter(fp1_particles,(n_decays,n_rows))
    E2_filter = make_filter(fp2_particles,(n_decays,n_rows))
    E3_filter = make_filter(fp3_particles,(n_decays,n_rows))
    E4_filter = make_filter(fp4_particles,(n_decays,n_rows))
    E5_filter = make_filter(fp5_particles,(n_decays,n_rows))
    # find the distance of each particle traveled in plate 1
    pop1_plate1_distance = np.sqrt(sum(fp1 - pos1)**2)
    pop2_plate1_distance = np.sqrt(sum(fp2 - pos2)**2)
    # find the distance of each particle traveled in plate 2
    pos4_P2 = fp2.copy()
    pos4_P2[fp4 == 0] = 0
    pop4_plate2_distance = np.sqrt(sum(fp4 - proj4)**2)
    pos5_P2 = fp4.copy()
    pos5_P2[fp5 == 0] = 0
    pop5_plate2_distance = np.sqrt(sum(fp5 - proj4)**2)
    # Find ED Plate 1 by pop1 and pop2
    ED1_P1 = -SEL.SpecificEnergyLoss(E1, material_values, pop1_plate1_distance)
    ED1_P1[E1_filter==0]=0
    ED2_P1 = -SEL.SpecificEnergyLoss(E2, material_values, pop2_plate1_distance)
    ED2_P1[E2_filter==0]=0
    # find RE of pop1 after passing through P1
    RE1_P1 = E1 - ED1_P1
    # Find RE of pop 3 and pop 4
    RE4_P1 = E2 - ED2_P1
    RE4_P1[E4_filter==0]=0
    RE3_P1 = E2 - ED2_P1
    RE3_P1[E3_filter==0]=0
    # Find ED in Plate 2 by pop4
    ED4_P2 = -SEL.SpecificEnergyLoss(E4, material_values, pop4_plate2_distance)
    ED5_P2 = -SEL.SpecificEnergyLoss(E5, material_values, pop5_plate2_distance)
    ED4_P2[E4_filter==0]=0
    ED5_P2[E5_filter==0]=0
    # Find RE of pop5
    RE5_P2 =  RE4_P1 - ED4_P2
    RE = [RE1_P1, RE3_P1, RE4_P1, RE5_P2]
    ED = [ED1_P1,ED2_P1,ED4_P2,ED5_P2]
    return RE, ED


def Conservation_of_Energy(E, RE,ED):
    RE1_P1, RE3_P1, RE4_P1, RE5_P2 = RE
    ED1_P1,ED2_P1,ED4_P2,ED5_P2 = ED
    # Energy conserved in Plate 1
    P1_Energy_Deposited = ED1_P1 + ED2_P1
    P2_Energy_Deposited = ED4_P2 + ED5_P2
    Not_Deposited = RE1_P1 + RE3_P1 + RE5_P2
    Total_Energy = P1_Energy_Deposited + P2_Energy_Deposited + Not_Deposited
    Energy_Cons = E - Total_Energy
    if sum(sum(Energy_Cons)) != 0:
        Titles = ['P1 Energy Deposited','P2 Energy Deposited','Not Deposited','Total Energy','Initial Energy','Energy Conservation']
        Energies = np.vstack([P1_Energy_Deposited,P2_Energy_Deposited,Not_Deposited,Total_Energy,E,Energy_Cons])
        Energies = Energies.reshape(len(Titles),E.size)
        print('Energy Not Conserved')
        write(Titles,Energies)
        print('csv file updated')
    return P1_Energy_Deposited,P2_Energy_Deposited,Not_Deposited


'''
PLOTS
'''
def vector(proja, projb,ax):
    xa, ya, za = remove_zeros(proja)
    xb, yb, zb = remove_zeros(projb)
    x = np.hstack([xa,xb]).reshape(2,len(xa))
    y = np.hstack([ya,yb]).reshape(2,len(xa))
    z = np.hstack([za,zb]).reshape(2,len(xa))
    # plate 1 subplots xz, xy, yz
    plt.suptitle(f'(Vacuum Vector, Electrons = {n_particles})')
    # xz 
    ax[0].plot(z,x,'--',color = 'black')
    # xy 
    ax[1].plot(x,y,'--',color = 'black')
    # yz 
    ax[2].plot(z,y,'--',color = 'black')

def plate1(ax):
    #xz
    ax[0].vlines(x= 0, ymin= -dim[1]/2, ymax = dim[1]/2, color = 'orange')
    ax[0].vlines(x= dim[2], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax[0].vlines(x= -dim[3], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax[0].hlines(y= -dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax[0].hlines(y= dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    #yx
    ax[1].vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[1].vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax[1].hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax[1].hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    #yz
    ax[2].vlines(x= 0, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'orange')
    ax[2].vlines(x= -dim[3], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[2].vlines(x= dim[2], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[2].hlines(y= -dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax[2].hlines(y= dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')

def plate2(ax):
    #xz
    ax[0].vlines(x= dim[2]+dim[5], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax[0].vlines(x= dim[2]+dim[5]+dim[4], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax[0].hlines(y= -dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
    ax[0].hlines(y= dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
    #yx
    ax[1].vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[1].vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax[1].hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax[1].hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    #yz
    ax[2].vlines(x= dim[2]+dim[5], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[2].vlines(x= dim[2]+dim[5]+dim[4], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax[2].hlines(y= -dim[0]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
    ax[2].hlines(y= dim[0]/2, xmin=dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')

def scatter(proj,axs,Lc,Lu,c):
    # xz
    axs[0].scatter(proj[2]*Lc,proj[0]*Lc,color = c)
    # yx
    axs[1].scatter(proj[0]*Lc,proj[1]*Lc,color = c)
    # yz
    axs[2].scatter(proj[2]*Lc,proj[1]*Lc,color = c)
    # labels
    axs[0].set_xlabel(f'z {Lu}')
    axs[1].set_xlabel(f'x {Lu}')
    axs[2].set_xlabel(f'z {Lu}')
    axs[0].set_ylabel(f'x {Lu}')
    axs[1].set_ylabel(f'y {Lu}')
    axs[2].set_ylabel(f'y {Lu}')

def HM_scatter(E,proj,axs,c,Ec,Lc,Lu):
    cm = mpl.cm.get_cmap(c)
    # xz
    axs[0].scatter(proj[2]*Lc,proj[0]*Lc,c = E*Ec, cmap = cm)
    # yx
    axs[1].scatter(proj[0]*Lc,proj[1]*Lc,c = E*Ec, cmap = cm)
    # yz
    sc = axs[2].scatter(proj[2]*Lc,proj[1]*Lc,c = E*Ec, cmap = cm)
    # labels
    axs[0].set_xlabel(f'z {Lu}')
    axs[1].set_xlabel(f'x {Lu}')
    axs[2].set_xlabel(f'z {Lu}')
    return sc

def E_hist(E,ED_P1,ED_P2,ND,Ec,Eu,n_particles):
    # energy information
    E = remove_zeros(E,n_decays*n_rows)*Ec
    ED_P1 = remove_zeros(ED_P1,n_decays*n_rows)*Ec
    ED_P2 = remove_zeros(ED_P2,n_decays*n_rows)*Ec
    ND = remove_zeros(ND,n_decays*n_rows)*Ec
    # only care about energy below 3keV
    E_3keV = E[E<=3]
    ED_P1_3keV = ED_P1[ED_P1<=3]
    ED_P2_3keV = ED_P2[ED_P2<=3]
    ND_3keV = ND[ND<=3]
    # find number of particles
    n_P1 = ED_P1.size
    n_P2 = ED_P2.size
    n_ND = ND.size
    # plots
    fig = plt.figure(figsize=(11,4))
    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)
    v = sum(E - np.mean(E)**2)/(len(E))
    # not normalized
    # ax1.hist(n1,bins1,color='red')
    # ax2.hist(n2,bins2,color='blue')
    # ax3.hist(n3,bins3,color='lime')
    # normalized
    # n1, b1 = np.histogram(ED_P1, bins = 100,range = (0,np.max(E)*Ec))
    # n2, b2 = np.histogram(ED_P2, bins = 10,range = (0,np.max(E)*Ec))
    # n3, b3 = np.histogram(ND, bins = 10,range = (0,np.max(E)*Ec))
    # bins1 = (b1[1:] + b1[:-1])/2
    # bins2 = (b2[1:] + b2[:-1])/2
    # bins3 = (b3[1:] + b3[:-1])/2
    ax1.hist(ED_P1_3keV,100,color='red')
    ax2.hist(ED_P2_3keV,100,color='blue')
    ax3.hist(ND_3keV,100,color='lime')
    ax4.hist(E_3keV,100,color='deeppink')
    ax1.set_xlabel(f'Energy {Eu}')
    ax2.set_xlabel(f'Energy {Eu}')
    ax3.set_xlabel(f'Energy {Eu}')
    ax4.set_xlabel(f'Energy {Eu}')
    ax1.set_ylabel('Intensity')
    ax1.set_title(f'Energy Detected P1 ({n_P1})')
    ax2.set_title(f'Energy Detected P2 ({n_P2})')
    ax3.set_title(f'Energy Not Detected ({n_ND})')
    ax4.set_title(f'Initial Energy ({n_particles})')
    plt.tight_layout()

'''
Position Graphs
'''
def Single_Plate(p,fsize,Lc,Lu,n_particles):
    p = np.vstack(p)
    p = p.reshape(5,3,n_decays,n_rows)
    p1 = remove_zeros(p[0],(n_decays,n_rows))
    p2 = remove_zeros(p[1],(n_decays,n_rows))
    p3 = remove_zeros(p[2],(n_decays,n_rows))
    p4 = remove_zeros(p[3],(n_decays,n_rows))
    p5 = remove_zeros(p[4],(n_decays,n_rows))
    # Plate 1
    fig1 = plt.figure(figsize=fsize)
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133) 
    axs = [ax1,ax2,ax3]
    plate1(axs)
    scatter(p1,axs,Lc,Lu,'blue')
    scatter(p2,axs,Lc,Lu,'red')
    ax1.set_ylabel('Plate1')
    ax2.set_title(f'Electrons = {n_particles}')
    # ax3.legend(['pop1','pop2'],loc = 'upper right',labelcolor = ['blue','red'])
    a1 = 0.01 #x
    a2 = 0.01 #y
    a3 = dim[1]/2*(dim[4]/dim[1]) #z
    ax1.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax1.set_xlim([-dim[3]-a3,dim[2]+a3])
    ax2.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax2.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_xlim([-dim[3]-a3,dim[2]+a3])
    plt.tight_layout()
    # Plate 2
    fig1 = plt.figure(figsize=fsize)
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133) 
    axs = [ax1,ax2,ax3]
    plate2(axs)
    scatter(p3,axs,Lc,Lu,'deeppink')
    scatter(p4,axs,Lc,Lu,'orange')
    scatter(p5,axs,Lc,Lu,'cyan')
    ax1.set_ylabel('Plate2')
    ax2.set_title(f'Electrons = {n_particles}')
    # ax3.legend(['pop3','pop4','pop5'],loc = 'upper right',labelcolor = ['deeppink','orange','cyan'])
    ax1.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax1.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    ax2.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax2.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    plt.tight_layout()

'''
Heat Maps
'''

def Single_Plate_HM(E_array,RE,ED,fps,initial_pos,proj,fsize,c,Ec,Lc,Eu,Lu,n_particles):
    # remove zeros from Energy arrays
    ED1_P1 = remove_zeros(ED[0],(n_decays*n_rows))
    ED2_P1 = remove_zeros(ED[1],(n_decays*n_rows))
    ED4_P2 = remove_zeros(ED[2],(n_decays*n_rows))
    ED5_P2 = remove_zeros(ED[3],(n_decays*n_rows)) # not on the graph
    RE1_P1 = remove_zeros(RE[0],(n_decays*n_rows))
    RE3_P1 = remove_zeros(RE[1],(n_decays*n_rows))
    RE4_P1 = remove_zeros(RE[2],(n_decays*n_rows)) # not on the graph
    RE5_P2 = remove_zeros(RE[3],(n_decays*n_rows))
    # remove zeros form position arrays
    fp1 = remove_zeros(fps[0],(n_decays*n_rows))
    fp2 = remove_zeros(fps[1],(n_decays*n_rows))
    fp3 = remove_zeros(fps[2],(n_decays*n_rows))
    fp4 = remove_zeros(fps[3],(n_decays*n_rows))
    fp5 = remove_zeros(fps[4],(n_decays*n_rows))
    # Initialize Plate 1 
    fig1 = plt.figure(figsize=fsize)
    ax1 = fig1.add_subplot(131)
    ax2 = fig1.add_subplot(132)
    ax3 = fig1.add_subplot(133) 
    axs = [ax1,ax2,ax3]
    plate1(axs)
    # scatter plots 
    HM_scatter(ED1_P1,fp1,axs,c,Ec,Lc,Lu)
    HM_scatter(RE1_P1,fp1,axs,c,Ec,Lc,Lu)
    sc2 = HM_scatter(ED2_P1,fp2,axs,c,Ec,Lc,Lu)
    # colorbar
    cbaxes = fig1.add_axes([0.1,0.03,0.8,0.01]) 
    cbar = plt.colorbar(sc2,cax = cbaxes, orientation = 'horizontal')
    cbar.set_label(f'Energy ({Eu})')
    # Labels
    ax1.set_ylabel('Plate1')
    ax2.set_title(f'Electrons = {n_particles}')
    a1 = 0.01 #x
    a2 = 0.01 #y
    a3 = dim[1]/2*(dim[4]/dim[1]) #z
    ax1.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax1.set_xlim([-dim[3]-a3,dim[2]+a3])
    ax2.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax2.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_xlim([-dim[3]-a3,dim[2]+a3])
    plt.tight_layout()
    # Initialize Plate 1 
    fig2 = plt.figure(figsize=fsize)
    ax1 = fig2.add_subplot(131)
    ax2 = fig2.add_subplot(132)
    ax3 = fig2.add_subplot(133) 
    axs = [ax1,ax2,ax3]
    plate2(axs)
    # scatter plots
    HM_scatter(RE3_P1,fp3,axs,c,Ec,Lc,Lu)
    sc4 = HM_scatter(ED4_P2,fp4,axs,c,Ec,Lc,Lu)
    HM_scatter(RE5_P2,fp5,axs,c,Ec,Lc,Lu)
    # colorbar
    cbaxes = fig1.add_axes([0.1,0.03,0.8,0.01]) 
    cbar = plt.colorbar(sc4,cax = cbaxes, orientation = 'horizontal')
    cbar.set_label(f'Energy ({Eu})')
    # Labels
    ax1.set_ylabel('Plate2')
    ax2.set_title(f'Electrons = {n_particles}')
    ax1.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax1.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    ax2.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax2.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    plt.tight_layout()

fps = Final_Positions(proj,E_array)
Conservation_of_Mass(E_array,fps)
initial_pos = Initial_Positions(proj)
RE,ED = Final_Energy(E_array,proj,fps, initial_pos)
P1_ED,P2_ED,ND = Conservation_of_Energy(E,RE,ED)
Single_Plate(fps,(9,3),1,'mm',n_particles) #p,fsize,Lc,Lu,n_particles
# Single_Plate_HM(E_array,RE,ED,fps,initial_pos,proj,(9,3),'plasma',1e-3,1,'keV','mm',n_particles)
# E_hist(E,P1_ED,P2_ED,ND,1e-3,'keV',n_particles)
plt.show()
print('END RUN')