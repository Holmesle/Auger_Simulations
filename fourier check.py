import matplotlib.pyplot as plt
import numpy as np


x = np.linspace(-np.pi,np.pi,int(1e5))
f = 0
for k in range(101):
    f += (2/np.pi)*np.sin(k*x)


plt.plot(x,f)
plt.show()
