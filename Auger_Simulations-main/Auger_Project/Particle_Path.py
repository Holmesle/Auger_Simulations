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
n_decays = 1000
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
source_to_inner_surface = 0.005e-3
source_to_outer_surface = 0.005e-3
thickness = source_to_inner_surface  + source_to_outer_surface
diameter = 0.4e-3  # 0.4 mm   # detector diameter (z)

'''
Required in functions
'''
dim, pos, angle, pos_matrices, angle_matrices = GM.Parallel_Plate(E, height, width, source_to_inner_surface, thickness, diameter)

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

def make_filter(v, shape):
    v= v.copy()
    if v.ndim > 2:
        v = sum(v)
    v[v==0]=0
    v[v!=0]=1
    v = v.reshape(shape)
    return v

'''
Determine the Projected Positons of all Populations 

for Reference:
dim = [height, width, innerplate_to_source,source_to_outerplate, thickness, diameter]
pos = [x0,y0,z0]
angle = [phi,theta]
'''
def projection(pos_matrices,angle_matrices,v_mag,n_decays):
    # source vector components
    x0 = pos_matrices[0:n_decays]
    y0 = pos_matrices[n_decays:n_decays*2]
    z0 = pos_matrices[n_decays*2:n_decays*3]
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
    proj1 = projection(pos_matrices, angle_matrices, -dim[3],n_decays)
    proj2 = projection(pos_matrices, angle_matrices, dim[2],n_decays)
    proj3 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]+dim[4]),n_decays)
    proj4 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]),n_decays)
    proj5 = projection(pos_matrices, angle_matrices, (dim[2]+dim[5]+dim[4]),n_decays)
    # sort out non-physical directions
    away = np.vstack([angle_matrices[n_decays:n_decays*2] > np.pi/2]*3)
    # use pop 3 because same distance as pop4 and nothing is filtered yet
    # pop5 doesnt need to be in plate after so long as it went through the plate
    outside_x = np.vstack([abs(proj3[0:n_decays]) < dim[1]/2]*3)
    outside_y = np.vstack([abs(proj3[n_decays:n_decays*2]) < dim[0]/2]*3)
    # apply angle conditions
    proj1[away] = 0          # pop1 goes away from P2
    proj2[~away] = 0        # pop2,3,4,5 go towards P2
    proj3[~away] = 0
    proj4[~away] = 0
    proj5[~away] = 0  
    # apply position conditions 
    proj4[~outside_x] = 0   # pop 3, 4 interact with P2
    proj4[~outside_y] = 0
    proj5[~outside_x] = 0   # pop 3, 4 interact with P2
    proj5[~outside_y] = 0
    proj3[proj4!=0] = 0
    # any zeros in E are empty place holders, carries through all pops
    no_particle = E == 0
    # no_particle = no_particle.reshape(n_decays*n_rows)
    no_particles  = np.vstack([no_particle]*3)
    # if energy is zero there is no particle there. 
    proj1[no_particles] = 0
    proj2[no_particles] = 0
    proj3[no_particles] = 0
    proj4[no_particles] = 0
    proj5[no_particles] = 0
    # conservation of matter
    # any particles with position values replace with ones and zeros
    p1_filter = make_filter(proj1,(3,n_decays*n_rows))
    p2_filter = make_filter(proj2,(3,n_decays*n_rows))
    p3_filter = make_filter(proj3,(3,n_decays*n_rows))
    p4_filter = make_filter(proj4,(3,n_decays*n_rows))
    p5_filter = make_filter(proj5,(3,n_decays*n_rows))
    # check conservation of matter
    detected_P1 = 2*sum(p1_filter)/3 + sum(p2_filter)/3
    detected_P2 = (sum(p4_filter)/3 + sum(p5_filter)/3)/2
    not_detected = sum(p3_filter)/3
    all_particles = (detected_P1 + detected_P2 + not_detected)/2
    n_rows = len(E[0,:])
    E_flat = E.reshape(n_rows*n_decays)
    E_flat_reduced = E_flat[E_flat!=0]
    n_particles = len(E_flat_reduced)
    if sum(all_particles) != n_particles:
        print('Error: mass not conserved')
        print('all:',all_particles)
        print('all sum:',sum(all_particles))
        print('n_particles: ',n_particles)
        print('n_decays: ',n_decays)
        print('n_rows: ',n_rows)
        # print('P1',detected_P1)
        # print('P2',detected_P2)
        # print('P3',not_detected)
        # print('proj3')
        # print(p3_filter)
    proj = np.vstack([proj1,proj2,proj3,proj4,proj5])
    # Energy associated with particle (no energy loss yet)
    E1 = E.copy()
    E2 = E.copy()
    E3 = E.copy()
    E4 = E.copy()
    E5 = E.copy()
    E1[sum(p1_filter).reshape(n_decays,n_rows) == 0] = 0
    E2[sum(p2_filter).reshape(n_decays,n_rows) == 0] = 0
    E3[sum(p3_filter).reshape(n_decays,n_rows) == 0] = 0
    E4[sum(p4_filter).reshape(n_decays,n_rows) == 0] = 0
    E5[sum(p5_filter).reshape(n_decays,n_rows) == 0] = 0
    E_array = np.vstack([E1,E2,E3,E4,E5])
    return proj, E_array


'''
Check Positions
'''
proj, E_array = Paths(E,n_decays,n_rows,pos_matrices,angle_matrices,dim)
proj = proj.reshape(5,3*n_decays,n_rows)

def xyz_scatter(E_array,proj,n_decays,n_rows,n_particles,dim):
    cm = mpl.cm.get_cmap('plasma')
    # unpack values
    proj1 = proj[0].reshape(3,n_decays,n_rows)
    proj2 = proj[1].reshape(3,n_decays,n_rows)
    proj3 = proj[2].reshape(3,n_decays,n_rows)
    proj4 = proj[3].reshape(3,n_decays,n_rows)
    proj5 = proj[4].reshape(3,n_decays,n_rows)
    E1 = E_array[0:n_decays]
    E2 = E_array[n_decays:n_decays*2]
    E3 = E_array[n_decays*2:n_decays*3]
    E4 = E_array[n_decays*3:n_decays*4]
    E5 = E_array[n_decays*4:n_decays*5]
    proj1 = remove_zeros(proj1,n_decays*n_rows)
    proj2 = remove_zeros(proj2,n_decays*n_rows)
    proj3 = remove_zeros(proj3,n_decays*n_rows)
    proj4 = remove_zeros(proj4,n_decays*n_rows)
    proj5 = remove_zeros(proj5,n_decays*n_rows)
    E1 = remove_zeros(E1,n_decays*n_rows)
    E2 = remove_zeros(E2,n_decays*n_rows)
    E3 = remove_zeros(E3,n_decays*n_rows)
    E4 = remove_zeros(E4,n_decays*n_rows)
    E5 = remove_zeros(E5,n_decays*n_rows)
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
    ax1.vlines(x= 0, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'orange')
    ax1.vlines(x= -dim[3], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.vlines(x= dim[2], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax1.hlines(y= -dim[1]/2, xmin= -dim[3], xmax = dim[2], color = 'k')
    ax1.hlines(y= dim[1]/2, xmin= -dim[3], xmax = dim[2], color = 'k')
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
    ax3.vlines(x= 0, ymin= -dim[0]/2, ymax = dim[0]/2, color = 'orange')
    ax3.vlines(x= -dim[3], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.vlines(x= dim[2], ymin= -dim[0]/2, ymax = dim[0]/2, color = 'k')
    ax3.hlines(y= -dim[0]/2, xmin= -dim[3], xmax = dim[2], color = 'k')
    ax3.hlines(y= dim[0]/2, xmin= -dim[3], xmax = dim[2], color = 'k')
    ax3.set_title('z vs y')
     # zx P2
    x345,y345,z345  = np.concatenate([proj5[0],proj4[0],proj3[0]]), np.concatenate([proj5[1],proj4[1],proj3[1]]), np.concatenate([proj5[2],proj4[2],proj3[2]])
    E345 = np.concatenate([E3,E4,E5])
    ax4.scatter(z345,x345,c = E345*Ec, cmap = cm)
    ax4.vlines(x= dim[2]+dim[5]+dim[4], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax4.vlines(x= dim[2]+dim[5], ymin= -dim[1]/2, ymax = dim[1]/2, color = 'k')
    ax4.hlines(y= -dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
    ax4.hlines(y= dim[1]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[5]+dim[4], color = 'k')
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
    ax6.hlines(y= -dim[0]/2, xmin= dim[2] + dim[5], xmax = dim[2]+dim[4]+dim[5], color = 'k')
    ax6.hlines(y= dim[0]/2, xmin= dim[2]+dim[5], xmax = dim[2]+dim[4]+dim[5], color = 'k')
    cbaxes = fig1.add_axes([0.1,0.03,0.8,0.01]) 
    cbar = plt.colorbar(sc2,cax = cbaxes, orientation = 'horizontal')
    cbar.set_label(f'Energy ({E_units})')
    a1 = 0.01 #x
    a2 = 0.01 #y
    a3 = dim[1]/2*(dim[4]/dim[1]) #z
    ax1.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax1.set_xlim([-dim[3]-a3,dim[2]+a3])
    ax2.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax2.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax3.set_xlim([-dim[3]-a3,dim[2]+a3])
    a1 = 0.01 #x
    a2 = 0.01 #y
    a3 = dim[1]/2*(dim[4]/dim[1]) #z
    ax4.set_ylim([-dim[1]/2-a1,dim[1]/2+a1])
    ax4.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    ax5.set_xlim([-dim[1]/2-a1,dim[1]/2+a1])
    ax5.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax6.set_ylim([-dim[0]/2-a2,dim[0]/2+a2])
    ax6.set_xlim([dim[5]+dim[2]-a3,dim[5]+dim[2]+dim[4]+a3])
    plt.tight_layout()

# xyz_scatter(E_array,proj,n_decays,n_rows,n_particles,dim)
# plt.show()

