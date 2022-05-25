import Parallel_Plate as PP

'''
This module requires dimension inputs from the detector module, uses the 
Track_Length and Final_Positions module to determine distance traveled 
through each plate and whether the particles hit the detector or remain 
undetected. In atempts to account for Secondary_Radiations or Interactions
I've made filler modules to be adapted later.

-> Geometry -> Secondary_Radiations -> Track_Length    -> Detector
            -> Interactions         -> Final Position

'''
height = 40e-3  # width (x) 40mm
width = 40e-3  # height (y)  40mm
source_to_inner_surface = 0.005e-3
source_to_outer_surface = 0.005e-3
thickness = source_to_inner_surface  + source_to_outer_surface
diameter = 0.4e-3  # 0.4 mm   # detector diameter (z)


def Parallel_Plate(E, height, width, source_to_inner_surface, source_to_outer_surface, diameter):
    thickness = source_to_inner_surface + source_to_outer_surface
    pos, angle, pos_matrices, angle_matrices = PP.Initial_Conditions(
        E, height, width)
    dim = [height, width, source_to_inner_surface,source_to_outer_surface,thickness,diameter]
    return dim, pos, angle, pos_matrices, angle_matrices


def Cylindrical(E, height, width, thickness, outerplate_to_source):
    pass
