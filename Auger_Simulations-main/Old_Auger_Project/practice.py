import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

# y = a1 + a2x + a3x**2 + ... + amx**m-1
x = np.arange(0,100,0.1)
n = len(x)
a = np.ones(3)
y = 3 + 2*x + x**2

def Normal_Eqns(a):
    ysum, xysum, x2ysum = sum(y), sum(x*y), sum(x**2*y)
    eqn1 = n*a[0] + a[1]*sum(x) + a[2]*sum(x**2)
    eqn2 = a[0]*sum(x) + a[1]*sum(x**2) + a[2]*sum(x**3)
    eqn3 = a[0]*sum(x**2) + a[1]*sum(x**3) + a[1]*sum(x**4)
    return[ysum - eqn1, xysum - eqn2, x2ysum - eqn3]

sol = optimize.root(Normal_Eqns, a, method = 'lm' )
a = sol.x
print(a)
y_curve = a[0] + a[1]*x + a[2]*x**2

plt.plot(x,y,color = 'black', label = 'True Eqn')
plt.plot(x, y_curve,'--', color ='red', label ='Curve Fitting')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.show()