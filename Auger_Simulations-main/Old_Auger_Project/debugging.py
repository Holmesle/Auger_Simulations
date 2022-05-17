import numpy as np
import matplotlib.pylab as plt
'''
# Wrong Shapes
'''

# a = np.array([1, 0, 0])
# b = np.array([0, 4, 0, 0])
# c = np.array([0, 0, 5])

# def Area_of_Parallelepiped(a,b,c):
#     a_cross_b = np.cross(a,b)
#     a_cross_b_dot_c = np.dot(a_cross_b,c)
#     return a_cross_b_dot_c

# print(Area_of_Parallelepiped(a, b, c))

''' 
Indexing error
'''
# there are two errors 

x = np.linspace(0,10,100)

def approx_derivative(x,f):
    # a likely familiar numerical approx
    # f'(x) ~ (f(x + h) - f(x))/h
    f_prime = np.zeros(len(f))
    for i in range(0,len(x)):
        h = x[i + 1] - x[i]
        f_prime[i] = (f[i+1] - f[i])/h

        



plt.figure()
plt.title('Linear')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.plot(x, f1,color = 'black',label = 'f(x)')
plt.plot(mid_x_f1, avg_y_f1, '*', color='lime', label='delta f(x)')
plt.legend()

plt.figure()
plt.title('Exponential')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.plot(x, f2, color='black', label='f(x)')
plt.plot(mid_x_f2, avg_y_f2, '*', color='red', label='delta f(x)')
plt.legend()
plt.show()

