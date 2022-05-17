import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import Geometry as GM
import Isotope_Source as IS
import Materials as MT

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

Pop1 = pop that is directed away from P2 on plate 1
Pop2 = pop that is directed toward P2 on P1
Pop3 = pop that is directed towards P2 but escapes
Pop4 = pop that is detected by P2
Pop5 = pop that passes through P2 (assuming high enough energy)

Conservation of Mass:
pop1 + pop 2 = pop3 + pop 4
(pop 5 is at some point a part of pop 4, more important in next modules)

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

'''
Required in functions
'''
dim, pos, angle, pos_matrices, angle_matrices = GM.Parallel_Plate(E, height, width, source_displacement, thickness, diameter)

'''
Determine the Projected Positons of all Populations 

for Reference:
dim = [height, width, innerplate_to_source,source_to_outerplate, thickness, diameter]
pos = [x0,y0,z0]
angle = [phi,theta]
'''
def projection(pos_matrices,angle_matrices,v_mag):
    # source vector components
    x0 = pos_matrices[0:2]
    y0 = pos_matrices[2:4]
    z0 = pos_matrices[4:6]
    # momentum directions
    theta = angle_matrices[0]
    phi = angle_matrices[1]
    # initiate vector component arrays
    x,y,z = np.zeros(np.shape(x0)),np.zeros(np.shape(x0)),np.zeros(np.shape(x0))
    # projection to length v_mag
    x += x0 + v_mag*np.cos(phi)*np.tan(theta)
    y += y0 + v_mag*np.sin(phi)*np.tan(theta)
    z += z0 + v_mag
    # concatneate
    r = np.vstack([[x][0],[y][0],[z][0]])
    return r

def Paths(E,n_decays,n_rows,pos_matrices,angle_matrices,dim):
    # project all electrons on each surface (maximal distances)
    proj_Pop1 = projection(pos_matrices, angle_matrices, -dim[3])
    proj_Pop2 = projection(pos_matrices, angle_matrices, dim[2])
    proj_Pop3 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]))
    proj_Pop4 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]))
    proj_Pop5 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]+dim[4]))
    # sort out non-physical directions
    up = np.vstack([angle_matrices[0:2] > np.pi/2]*3)
    inside_x = np.vstack([abs(proj_Pop3[0:2]) < dim[1]/2]*3)
    inside_y = np.vstack([abs(proj_Pop3[2:4]) < dim[0]/2]*3)
    # apply angle conditions
    proj_Pop1[up] = 0          # pop1 goes away from P2
    proj_Pop2[~up] = 0        # pop2,3,4,5 go towards P2
    proj_Pop3[~up] = 0
    proj_Pop4[~up] = 0
    proj_Pop5[~up] = 0   
    # apply position conditions
    proj_Pop3[inside_x] = 0    # pop 3 misses P2
    proj_Pop3[inside_y] = 0   
    proj_Pop4[~inside_x] = 0   # pop 3, 4 interact with P2
    proj_Pop4[~inside_y] = 0
    proj_Pop5[~inside_x] = 0   # pop 3, 4 interact with P2
    proj_Pop5[~inside_y] = 0
    # any zeros in E are empty place holders, carries through all pops
    no_particle = E == 0
    # no_particle = no_particle.reshape(n_decays*n_rows)
    no_particles  = np.vstack([no_particle,no_particle,no_particle])
    # if energy is zero there is no particle there. 
    proj_Pop1[no_particles] = 0
    proj_Pop2[no_particles] = 0
    proj_Pop3[no_particles] = 0
    proj_Pop4[no_particles] = 0
    proj_Pop5[no_particles] = 0
    # conservation of matter
    pop1 = sum(proj_Pop1.reshape(3,n_decays*n_rows))
    pop2 = sum(proj_Pop2.reshape(3,n_decays*n_rows))
    pop3 = sum(proj_Pop3.reshape(3,n_decays*n_rows))
    pop4 = sum(proj_Pop4.reshape(3,n_decays*n_rows))
    pop5 = sum(proj_Pop5.reshape(3,n_decays*n_rows))
    # any particles with position values replace with ones and zeros
    pop1[pop1 != 0] = 1
    pop2[pop2 != 0] = 1
    pop3[pop3 != 0] = 1
    pop4[pop4 != 0] = 1
    pop5[pop5 != 0] = 1
    # check conservation of matter
    Mcons = pop1 + pop3 + pop4 + pop2 - pop5
    Mcons = Mcons.reshape(n_decays,n_rows)
    Mcons = Mcons[E!=0]
    if Mcons.all() != np.ones(np.shape(pop1.reshape(n_decays,n_rows)[E!=0])).all():
        print(Mcons)
        print(np.ones(np.shape(pop1)))
        print('Mass is not conserved')
    proj = np.vstack([proj_Pop1,proj_Pop2,proj_Pop3,proj_Pop4,proj_Pop5])
    # Energy associated with particle (no energy loss yet)
    E1 = E.copy()
    E2 = E.copy()
    E3 = E.copy()
    E4 = E.copy()
    E5 = E.copy()
    E1[pop1.reshape(n_decays,n_rows) == 0] = 0
    E2[pop2.reshape(n_decays,n_rows) == 0] = 0
    E3[pop3.reshape(n_decays,n_rows) == 0] = 0
    E4[pop4.reshape(n_decays,n_rows) == 0] = 0
    E5[pop5.reshape(n_decays,n_rows) == 0] = 0
    E_array = np.vstack([E1,E2,E3,E4,E5])
    return proj, E_array


'''
Check Positions
'''
proj, E_array = Paths(E,n_decays,n_rows,pos_matrices,angle_matrices,dim)
proj = proj.reshape(5,3*n_decays,n_rows)
def xyz_scatter(E_array,proj,n_decays,n_rows,n_particles,dim):
    cm = mpl.cm.get_cmap('plasma')
    print('proj = ', np.shape(proj))
    # unpack values
    proj1 = proj[0].reshape(3,n_decays,n_rows)
    proj2 = proj[1].reshape(3,n_decays,n_rows)
    proj3 = proj[2].reshape(3,n_decays,n_rows)
    proj4 = proj[3].reshape(3,n_decays,n_rows)
    proj5 = proj[4].reshape(3,n_decays,n_rows)
    E1 = E_array[0:2]
    E2 = E_array[2:4]
    E3 = E_array[4:6]
    E4 = E_array[6:8]
    E5 = E_array[8:10]
    # how many particles
    np1 = (len(E1[E1!=0]))
    np2 = (len(E2[E2!=0]))
    np3 = (len(E3[E3!=0]))
    np4 = (len(E4[E4!=0]))
    np5 = (len(E5[E5!=0]))
    proj1 = proj1[proj1!=0].reshape(3,np1)
    proj2 = proj2[proj2!=0].reshape(3,np2)
    proj3 = proj3[proj3!=0].reshape(3,np3)
    proj4 = proj4[proj4!=0].reshape(3,np4)
    proj5 = proj5[proj5!=0].reshape(3,np5)
    E1 = E1[E1!=0].reshape(np1)
    E2 = E2[E2!=0].reshape(np2)
    E3 = E3[E3!=0].reshape(np3)
    E4 = E4[E4!=0].reshape(np4)
    E5 = E5[E5!=0].reshape(np5)
    # plate 1 subplots xz,xy,yz
    fig1 = plt.figure(figsize=(9,6))
    ax1 = fig1.add_subplot(231)
    ax2 = fig1.add_subplot(232)
    ax3 = fig1.add_subplot(233)
    ax4 = fig1.add_subplot(234)
    ax5 = fig1.add_subplot(235)
    ax6 = fig1.add_subplot(236)
    plt.suptitle(f'(Electrons = {n_particles})')
    # xz P1
    Ec = 1e-3
    E_units = 'keV'
    z12,x12,y12  = np.concatenate([proj1[2],proj2[2]]), np.concatenate([proj1[0],proj2[0]]),np.concatenate([proj1[1],proj2[1]])
    E12 = np.concatenate([E1,E2])
    ax1.scatter(z12,x12,c = E12, cmap = cm)
    ax1.vlines(x= -dim[3], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.vlines(x= dim[2], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.hlines(y= -dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.hlines(y= dim[1]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax1.set_title('x vs z')
    ax1.set_ylabel('Plate1')
    # xy P1
    ax2.scatter(x12,y12,c = E12*Ec, cmap = cm)
    ax2.vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax2.vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax2.hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax2.set_title('y vs x')
    # zy P1
    ax3.scatter(z12,y12,c = E12*Ec, cmap = cm)
    ax3.vlines(x= -dim[3], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.vlines(x= dim[2], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.hlines(y= -dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.hlines(y= dim[0]/2, xmin= dim[2], xmax = -dim[3], color = 'k')
    ax3.set_title('z vs y')
     # zx P2
    x345,y345,z345  = np.concatenate([proj5[0],proj4[0],proj3[0]]), np.concatenate([proj5[1],proj4[1],proj3[1]]), np.concatenate([proj5[2],proj4[2],proj3[2]])
    E345 = np.concatenate([E3,E4,E5])
    ax4.scatter(z345,x345,c = E345*Ec, cmap = cm)
    ax4.vlines(x= dim[2]+dim[5]+dim[4], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax4.vlines(x= dim[2]+dim[5], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax4.hlines(y= -dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
    ax4.hlines(y= dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[3]+dim[5]+dim[4], color = 'k')
    ax4.set_ylabel('Plate2')
    # # xy P2
    ax5.scatter(x345,y345,c = E345*Ec, cmap = cm)
    ax5.vlines(x= dim[1]/2, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax5.vlines(x= -dim[1]/2, ymin= -dim[0]/2, ymax = dim[1]/2, color = 'k')
    ax5.hlines(y= -dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax5.hlines(y= dim[0]/2, xmin= -dim[1]/2, xmax = dim[1]/2, color = 'k')
    ax5.set_xlabel('    ')
    # # zy P2
    sc2 = ax6.scatter(z345,y345,c = E345*Ec, cmap = cm)
    ax6.vlines(x= dim[2]+dim[5]+dim[4], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax6.vlines(x= dim[2]+dim[5], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax6.hlines(y= -dim[0]/2, xmin= dim[2]+dim[5], xmax = dim[4]+dim[2]+dim[5], color = 'k')
    ax6.hlines(y= dim[0]/2, xmin= dim[2]+dim[5], xmax = dim[4]+dim[2]+dim[5], color = 'k')
    cbaxes = fig1.add_axes([0.1,0.03,0.8,0.01]) 
    cbar = plt.colorbar(sc2,cax = cbaxes, orientation = 'horizontal')
    cbar.set_label(f'Energy ({E_units})')
    plt.tight_layout()

xyz_scatter(E_array,proj,n_decays,n_rows,n_particles,dim)
plt.show()

