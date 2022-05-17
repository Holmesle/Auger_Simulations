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

def Parallel_Plate(E, height, width, innerplate_to_source, thickness, diameter):
    source_to_outerplate = abs(thickness - innerplate_to_source)
    pos, angle, pos_matrices, angle_matrices = PP.Initial_Conditions(
        E, height, width)
    dim = [height, width, innerplate_to_source,
           source_to_outerplate, thickness, diameter]
    return dim, pos, angle, pos_matrices, angle_matrices


def Cylindrical(E, height, width, thickness, outerplate_to_source):
    pass
