import numpy as np
import matplotlib.pyplot as plt

m = 1
h = 1
dr = 0.1
k = 1
l = 0
L = 10
E = 10

# Example 1

U = np.zeros(int(L/dr))
U[-2] = 1
V = np.zeros(int(L/dr))

for r in range(int(L/dr)):
    U[r-1] = 2*U[r] - U[r+1] + dr**2*(l*(l+1)/(r**2)+2*m*(V - E))*U[r]

# Example 2

U = np.zeros(int(L/dr))
U[-2] = 1

for r in range(int(L/dr)):
    V = k*r**2

for i in range(int(L/dr)):
    U[i-1] = 2*U[i] - U[i+1] + dr**2*(l*(l+1)/(r[i]**2)+2*m*(V[i] - E))*U[i]

# Example 3

U = np.zeros(int(L/dr))
U[-2] = 1
V = np.zeros(int(L/dr))
r = np.zeros(int(L/dr))

V = k*r**2

for i in range(int(L/dr)):
    U[i-1] = 2*U[i] - U[i+1] + dr**2*(l*(l+1)/(r[i]**2)+2*m*(V[i] - E))*U[i]


r = np.arange(0,L,dr)
V = k*r**2