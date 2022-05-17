import matplotlib.pylab as plt
import numpy as np

plt.figure()
plt.xlim([-10,10])
plt.ylim([-10, 10])

x = np.arange(-5,5,0.01)
r = 6
print(r**2 - x**3)
semi_circle = np.sqrt(r**2 - x**2)
trunk_r = -2*(x-3)**2 -1
plt.plot(x,semi_circle,color = 'black')
plt.plot(x, trunk_r, color='black')
plt.show()
