'''
RL circuit data
'''
import numpy as np
import matplotlib.pyplot as plt

print('start')

with open('RL-data.txt','r') as file:
    line = file.readline()
    print(line)
file.close()

print('end')