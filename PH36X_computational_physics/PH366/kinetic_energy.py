from cProfile import label
import re
import numpy as np
import matplotlib.pyplot as plt

#### variables

m = 1                  # 9.109e-31
hbar = 1                  # 1.05e-34

L = 5
dx = 0.01
x = np.arange(0, L, dx)

N = len(x)                   # dim of matrix T

k = 1                     # spring const

'''
2. create a matrix that represents kinetic energy
'''


# T operator * psi function
def T(N):
    k = hbar**2/(2*m*dx**2)
    M = np.zeros([N,N])
    a0 =np.arange(N-1)
    a1 = np.arange(0,N)
    a2 = np.arange(1,N)
    M[a0, a0+1] = -1
    M[a1, a1] = 2
    M[a2,a2-1] = -1
    return k*M

# print(T(4))

'''
3. find eigenv values and vectors of T*psi matrix
'''

eigT_value, eigT_vectors = np.linalg.eigh(T(N))

'''
4. Plot 4 eigenfunctions (lowest energies)
'''
colors = ['midnightblue','blue','blueviolet','darkviolet','mediumorchid',
'violet','hotpink','deeppink','red','orangered','tomato','salmon',
'coral','lightsalmon','moccasin','gold','khaki','yellowgreen']

# fig, axs = plt.subplots(nrows=6, ncols=1)

# fig.set_size_inches(8,6)
# axs = axs.flatten()
# for i in range(6):
#     axs[i].plot(x, eigT_vectors[:,i], color=colors[i])
#     axs[i].set_ylabel(f'T*Ψ({i+1}*dx)')
#     axs[i].set_xlabel('x')
#     plt.tight_layout()


'''
5. Plot eigenvalues
'''
# plt.figure(figsize=(3,6))
# for i in range(6):
#     plt.plot([0, 1],[eigT_value[i], eigT_value[i]],color = colors[i],label = f'n = {i+1}')
#     plt.legend(loc = 'lower right')
# plt.show()

'''
Extra Fun
'''


def X(N,x):
    M = np.zeros([N, N])
    din = np.diag_indices_from(M)
    M[din] = x
    return M

print(X(5,np.linspace(0,L+1,5)))

V = k*X(N,x)**2*0.5

eigV_value, eigV_vectors = np.linalg.eigh(V)

fig, axs = plt.subplots(nrows=5, ncols=1)

fig.set_size_inches(6, 6)

for i in range(5):
    axs[i].plot(x, eigV_vectors[:, i], color=colors[i])
    axs[i].set_ylabel(f'Ψ({i+1}*dx)')
    axs[i].set_xlabel('x')
    axs[i].set_xlim([0, 0.07])
    plt.tight_layout()

plt.figure(figsize=(3, 6))
for i in range(4+1):
    plt.plot([0, 1], [eigV_value[i], eigV_value[i]],
             color=colors[i], label=f'n = {i+1}')
    plt.legend(loc='lower right')

# plt.show()

'''
Springy Fun
'''

H = T(N) + V

eigH_value, eigH_vectors = np.linalg.eigh(H)

fig, axs = plt.subplots(nrows= 5, ncols=1)

fig.set_size_inches(6,6)
axs = axs.flatten()

for i in range(5):
    axs[i].plot(x, eigH_vectors[:,i], color=colors[i])
    axs[i].set_ylabel(f'Ψ({i+1}*dx)')
    axs[i].set_xlabel('x')
    plt.tight_layout()

plt.figure(figsize=(3,6))
plt.plot(x, k*x**2*0.5,color = 'black')
for i in range(4+1):
    plt.plot([0, L],[eigH_value[i], eigH_value[i]],color = colors[i],label = f'n = {i+1}')
    plt.legend(loc = 'lower right')
plt.show()
