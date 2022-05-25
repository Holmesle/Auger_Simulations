from typing import Final
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Geometry as GM
import Isotope_Source as IS
import Materials as MT
import Particle_Path as PTH
import Specific_Energy_Loss as SEL
import Parallel_Plate as PP

'''
This module requires initial energy inputs from Isotope_Source and material 
properties form Materials module and outputs the final and deposited 
energies, and track lengths to the Detector module.

Isotope_Source -> Particle_Path -> Energy_deposited -> Detector
     Materials ->                  
      Geometry ->  

Variable Names Explained:
PopA = populated and detected in Plate 1 (P1)
PopB = population detected in Plate 2 (P2)
PopC = Population not detected by P2

Pop1 = pop that is directed away from P2 on P1
Pop2 = pop that is directed toward P2 on P1
Pop3 = pop that is directed towards P2 but escapes
Pop4 = pop that is detected by P2
Pop5 = pop that passes through P2 (assuming high enough energy)

Pop A = pop1 + pop2
Pop B = pop4 - pop5
Pop C = pop3 + pop 5

Conservation of Energy:
E - E(pop1) -E(pop2) - E3(pop3) - E4(pop4)

'''

'''
 Varaibles Declared in the Detector Module
'''
# Energy and Particles
n_decays = 2
isotope = 'I125'
m, E = IS.JANIS_properties(n_decays, isotope)
n_rows = len(E[0,:])
E_flat = E.reshape(n_rows*n_decays)
E_flat_reduced = E_flat[E_flat!=0]
n_particles = len(E_flat_reduced)
# Materials
material = 'Si'
material_values = MT.material_dict(material)  # material values
# Plate dimensions
height = 40e-3  # width (x) 40mm
width = 40e-3  # height (y)  40mm
source_displacement = 0.005e-3
diameter = 0.4e-3  # 0.4 mm   # detector diameter (z)
thickness = source_displacement*2  # thickness of the plate 0.02mm
distance_source_to_inner_plate2 = diameter + source_displacement
# Geometry
dim, pos, angle, pos_matrices, angle_matrices = GM.Parallel_Plate(E, height, width, source_displacement, thickness, diameter)
print('new run')

'''
Projection Vector Length

SEL requires s, so we need to determine the maximum possible 
traveling distance through the plates
'''
proj, E_array = PTH.Paths(E,n_decays,n_rows,pos_matrices,angle_matrices,dim)
proj = proj.reshape(5,3*n_decays,n_rows) # shape so that we can call each pop
proj1 = proj[0].reshape(3,n_decays,n_rows)
proj2 = proj[1].reshape(3,n_decays,n_rows)
proj3 = proj[2].reshape(3,n_decays,n_rows)
proj4 = proj[3].reshape(3,n_decays,n_rows)
proj5 = proj[4].reshape(3,n_decays,n_rows)
# magnitude of projected vector
mag1 = np.sqrt(proj1[0]**2 + proj1[1]**2 + proj1[2]**2)
mag2 = np.sqrt(proj2[0]**2 + proj2[1]**2 + proj2[2]**2)
mag3 = np.sqrt(proj3[0]**2 + proj3[1]**2 + proj3[2]**2)
mag4 = np.sqrt(proj4[0]**2 + proj4[1]**2 + proj4[2]**2)
mag5 = np.sqrt(proj5[0]**2 + proj5[1]**2 + proj5[2]**2)
E1 = E_array[0:2]
E2 = E_array[2:4]
E3 = E_array[4:6]
E4 = E_array[6:8]
E5 = E_array[8:10]
# find vectors from P1 to P2 within vacuum
vac = PTH.projection(proj2.reshape(3*n_decays,n_rows),angle_matrices,dim[5])
vac = vac.reshape(3,n_decays,n_rows)
vac1 = vac.copy()
vac2 = vac.copy()
vac3 = vac.copy()
vac4 = vac.copy()
vac5 = vac.copy()
vac1[proj1 ==0] = 0
vac2[proj2 ==0] = 0
vac3[proj3 ==0] = 0
vac4[proj4 ==0] = 0
vac5[proj5 ==0] = 0
# find distance traveled in vacuum (no energy is lost in this region)
vmag1 = np.sqrt(vac1[0]**2 + vac1[1]**2 + vac1[2]**2)
vmag2 = np.sqrt(vac2[0]**2 + vac2[1]**2 + vac2[2]**2)
vmag3 = np.sqrt(vac3[0]**2 + vac3[1]**2 + vac3[2]**2)
vmag4 = np.sqrt(vac4[0]**2 + vac4[1]**2 + vac4[2]**2)
vmag5 = np.sqrt(vac5[0]**2 + vac5[1]**2 + vac5[2]**2)
'''
Allowable Vector Length

using max_s and the enegies find the actual distance traveled
and the find the particles final positions
'''
# find max allowable distance based on inital energies
s1 = SEL.max_s(E1,n_rows,n_decays)*1e-7 #convert to milimeters
s2 = SEL.max_s(E2,n_rows,n_decays)*1e-7 
s3 = SEL.max_s(E3,n_rows,n_decays)*1e-7 
s4 = SEL.max_s(E4,n_rows,n_decays)*1e-7 
s5 = SEL.max_s(E5,n_rows,n_decays)*1e-7 
# all particles that do not meet conditions of pop are removed
s1[mag1 == 0] = 0
s2[mag2 == 0] = 0
s3[mag3 == 0] = 0
s4[mag4 == 0] = 0
s5[mag5 == 0] = 0

'''
Determine Actual Position Vectors
'''
# find remaining distances after plate 1
dif1 = s1 - mag1
dif2 = s2 - mag2
# find remaining distances after plate 2
dif3 = s3 - (mag3 - vmag3)
dif4 = s4 - (mag4 - vmag4)
dif5 = s5 - (mag5 - vmag5)
# if values are negative then the paticle stops in the plate
negP1a = dif1 < 0
negP1b = dif2 < 0
negP2 = dif5 < 0 
# need to adjust final position vectors in proj 2 and proj 4
stopsP1a = -s1.copy()
stopsP1b = -s2.copy()
stopsP2 = -s4.copy()
stopsP1a[~negP1a] = 0
stopsP1b[~negP1b] = 0
stopsP2[~negP2] = 0
# find the vectors of these magnitudes.
proj_adj1 = PTH.projection(pos_matrices,angle_matrices,stopsP1a).reshape(3,n_decays,n_rows)
proj_adj2 = PTH.projection(pos_matrices,angle_matrices,stopsP1b).reshape(3,n_decays,n_rows)
proj_adj4 = PTH.projection(pos_matrices,angle_matrices,stopsP2).reshape(3,n_decays,n_rows)
# replace particles that stop in P1 with s values
fp1 = proj1.copy()
fp2 = proj2.copy()
fp1[proj_adj1 != 0] = proj_adj1[proj_adj1 != 0 ]
fp2[proj_adj2 != 0] = proj_adj2[proj_adj2 != 0 ]
#remove empty place holders
fp1[proj1==0] = 0
fp2[proj2==0] = 0
print('dif1:',dif1[dif1!=0])
fp1z = fp1[2]
projz = proj1[0]
print('projz1:')
print(projz[projz!=0])
print(np.shape(proj1[0]))
print('fp1z:')
print(fp1z[fp1z!=0])
print(np.shape(fp1z[0]))
# final magnitudes
fm1 = np.sqrt(fp1[0]**2 + fp1[1]**2 + fp1[2]**2)
fm2 = np.sqrt(fp2[0]**2 + fp2[1]**2 + fp2[2]**2)
#Using the final magnitude determine the ED = -SEL
ED1 = -SEL.SpecificEnergyLoss(E1, material_values, fm1)
ED2 = -SEL.SpecificEnergyLoss(E2, material_values, fm2)
RE1 = E1 - ED1 
RE2 = E2 - ED2 

'''
Plot
'''
def xyz_cm(E_array,proj1,n_particles,dim,fig):
    cm = mpl.cm.get_cmap('plasma')
    # unpack values
    E1 = E_array[0:2]
    # how many particles
    np1 = (len(E1[E1!=0]))
    proj1 = proj1[proj1!=0].reshape(3,np1)
    E1 = E1[E1!=0].reshape(np1)
    # plate 1 subplots xz,xy,yz
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    plt.suptitle(f'(Electrons = {n_particles})')
    # xz P1
    Ec = 1e-3
    E_units = 'keV'
    z1,x1,y1  = proj1[2],proj1[0],proj1[1]
    ax1.scatter(z1,x1,c = E1*Ec, cmap = cm)
    ax1.vlines(x= -dim[3], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.vlines(x= dim[2], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.hlines(y= -dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.hlines(y= dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.set_title('x vs z')
    ax1.set_ylabel('Plate1')
    # xy P1
    ax2.scatter(x1,y1,c = E1*Ec, cmap = cm)
    ax2.vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax2.vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax2.hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.set_title('y vs x')
    # zy P1
    sc = ax3.scatter(z1,y1,c = E1*Ec, cmap = cm)
    ax3.vlines(x= -dim[3], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.vlines(x= dim[2], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.hlines(y= -dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.hlines(y= dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.set_title('z vs y')
    cbaxes = fig.add_axes([0.1,0.03,0.8,0.01]) 
    cbar = plt.colorbar(sc,cax = cbaxes, orientation = 'horizontal')
    cbar.set_label(f'Energy ({E_units})')
    plt.tight_layout()

def xyz_c(E_array,proj1,n_particles,dim,fig,c):
    # unpack values
    E1 = E_array[0:2]
    # how many particles
    np1 = (len(E1[E1!=0]))
    proj1 = proj1[proj1!=0].reshape(3,np1)
    E1 = E1[E1!=0].reshape(np1)
    # plate 1 subplots xz,xy,yz
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    print(np.shape(proj1))
    plt.suptitle(f'(Electrons = {n_particles})')
    # xz P1
    E_units = 'keV'
    z1,x1,y1  = proj1[2],proj1[0],proj1[1]
    ax1.scatter(z1,x1,color = c)
    ax1.vlines(x= -dim[3], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.vlines(x= dim[2], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.hlines(y= -dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.hlines(y= dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.set_title('x vs z')
    ax1.set_ylabel('Plate1')
    # xy P1
    ax2.scatter(x1,y1,color = c)
    ax2.vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax2.vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax2.hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.set_title('y vs x')
    # zy P1
    sc = ax3.scatter(z1,y1,color = c)
    ax3.vlines(x= -dim[3], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.vlines(x= dim[2], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.hlines(y= -dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.hlines(y= dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.set_title('z vs y')
    plt.tight_layout()

print(PP.Initial_Conditions(E,height,width)[0])
print(proj1)
# fig1 = plt.figure(figsize=(9,3))
# fig2 = plt.figure(figsize=(9,3))
# fig3 = plt.figure(figsize=(9,3))
# fig4 = plt.figure(figsize=(9,3))
# xyz_cm(ED1,fp1,n_particles,dim,fig1)
# xyz_cm(ED2,fp2,n_particles,dim,fig2)
# xyz_c(E1,proj1,n_particles,dim,fig3,'blue')
# plt.suptitle('Plate 1 Pop1')
# xyz_c(E2,proj2,n_particles,dim,fig4,'lime')
# plt.suptitle('Plate 1 Pop2')
# projx = proj1[0]
# print('dif1:',dif1[dif1!=0])
# print('projx1:',projx[projx!=0],np.shape(proj1[0]))
# print('fm1:',fm1[fm1!=0],np.shape(fm1))
# print('mag1:',mag1[mag1!=0],np.shape(mag1))
# print('s1:',s1[s1!=0],np.shape(s1))
# fp1x = fp1[0]
# fp1z = fp1[2]
# projz = proj1[0]
# print('dif1:',dif1[dif1!=0])
# print('projz1:',projz[projz!=0],np.shape(proj1[0]))
# print('fp1x:',fp1x[fp1x!=0])
# print('fp1z:',fp1z[fp1z!=0])
# plt.show()

'''
Energy Remaining

Use the energy deposited to determine if there is any remaining 
energy in population 5.
'''

'''
Check Conservation of Energy
'''


'''
Scatter Plots
'''

'''
Histograms
'''



