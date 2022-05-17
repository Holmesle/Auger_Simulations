import matplotlib.pylab as plt
import numpy as np
import pandas as pd

def material_dict(Element_Formula):
    data_filename = r"C:\Users\holme\OneDrive\Desktop\python_code\Auger_Simulations-main\Materials.csv"
    with open(data_filename, 'r') as myfile:
        d = dict(pd.read_csv(myfile, delimiter=',', skiprows=[0], index_col=0))
    return d[Element_Formula]

#### check ####
# k,J,Z,A,p = material_dict('Si')
# print(k)

