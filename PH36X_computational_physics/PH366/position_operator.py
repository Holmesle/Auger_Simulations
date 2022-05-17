import numpy as np
import matplotlib.pyplot as plt

'''
1. write a function that returns xnm
'''
L = 1
N = 10
dx = 0.01
x = np.arange(0,L,dx)

def phi(n,x):
    return np.sqrt(2/L)*np.sin(n*x*np.pi/L)

def Xnm(phi_n,phi_m,x):
    I = 0
    np.shape(phi_n)
    np.shape(phi_m)
    np.shape(x)
    for i in range(0,len(x)):
        I += phi_n[i]*x[i]*phi_m[i]*dx
    return I

'''
2. Create matrix for position operator
'''
def matrix(N):
    X = np.zeros([N,N])
    for i in range(1,N+1):
        for j in range(1,N+1):
            X[i-1,j-1] = Xnm(phi(i,x),phi(j,x),x)
    return X

X = matrix(N)

'''
3. Plot heat map of the position operator
'''

# plt.figure()
# c = plt.imshow(X)
# plt.colorbar(c)
# plt.xlabel('m')
# plt.ylabel('n')
# plt.show()

'''
1. solve for the eigenvalues and vectors of the position matrix
'''

eignevalues, eigenvectors = np.linalg.eig(X)


def vi(eigvector,phi_n):
    return sum(i * phi_n for i in eigvector)

colors = ['red','lime','blue']

# plt.figure()
# for i in range(1,4):
#     plt.plot(x,vi(eigenvectors[i-1],phi(i,x)),color = colors[i-1], label = f'eig vector {i}')
#     plt.plot([eignevalues[i-1], eignevalues[i-1]], [-3, 3],'--',color = colors[i-1], label = f'eig value {i}')
#     plt.legend()
# plt.show()
