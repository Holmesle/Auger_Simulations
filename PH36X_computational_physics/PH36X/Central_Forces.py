from tracemalloc import stop
import numpy as np
import matplotlib.pyplot as plt

'''
Variables & Parameters
'''

R = 1 # max radius
l = 0.5 # angular momentum quantum number
h = 1   # reduced plank const
m = 1   # mass
k = 0.5   # string const
E = 0.7071

'''
Arrays
'''
dr = 0.01
r = np.arange(0,R+0.01,dr)
# r = r[::-1]     #reverse order of r
u = np.zeros(r.size)
u[-2] = 5


'''
Guess
'''

E0 = (0+0.5)*h*np.sqrt(k/m)
E1 = (1+0.5)*h*np.sqrt(k/m)
E2 = (2+0.5)*h*np.sqrt(k/m)
E3 = (3+0.5)*h*np.sqrt(k/m)

print('E:',E0,E1,E2,E3)

'''
Potential Function
'''
def Vs(r):
    return 0.5*k*r**2

'''
Radial Wave Function
'''
# # For Loop
for i in range(r.size-2,-1,-1):
    # print(u[i-1])
    # print(u[i])
    # print(u[i+1])
    V = Vs(r[i])
    u[i-1] = 2*u[i] - u[i+1] + dr**2*(l*(l+1)/(r[i]**2) 
    + (2*m/(h**2)*(V-E)))*u[i]

# Slicing
# u1 = np.zeros(r.size)
# u1[2] = 2
# r1 = np.arange(0,R+0.01,dr)
# r1 = r[::-1]     #reverse order of r

# u1[2:] = 2*u1[1:-1] - u1[:-2] + dr**2*(l*(l+1)/(r1[1:-1]**2) 
# + (2*m/(h**2)*(V-np.ones(V.size)*E)))*u1[1:-1]
    
plt.plot(r,u,'-',color = 'red')
plt.plot(r,np.zeros(r.size),color = 'black')
# plt.plot(r1,u1,'--',color = 'blue')
plt.show()



