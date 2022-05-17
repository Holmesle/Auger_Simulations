from sys import builtin_module_names
import numpy as np
import matplotlib.pyplot as plt

b_num = 10
t_max = 100
dt = 0.1
t_range = np.arange(0,t_max,dt)
T = 1
m = 1
dx = 0.01


y = np.zeros(b_num,t_max)

for b in range(b_num):
    for t in range(len(t)):
        y[b,t] = 2*y[b,t] - y[b,t] + (T/m*dx)