'''
find mean x^2 for random 1 link w/ length L
'''
import numpy as np
import matplotlib.pyplot as plt

L = 1
N = 20  # number of angles
colors = ['red', 'gold', 'orange', 'lime', 'cyan',
          'blue', 'violet', 'pink', 'black', 'red']

def Links(L,N):
    thetas = (np.random.rand(N))*2*np.pi
    x_array = np.zeros(N)
    y_array = np.zeros(N)
    for i in range(1,N):
        new_x = L*np.cos(thetas[i])
        new_y = L*np.sin(thetas[i])
        old_x = x_array[i-1]
        old_y = y_array[i-1]
        x_array[i] = old_x + new_x
        y_array[i] = old_y+ new_y
    return x_array,y_array

x, y = Links(L,N)

mean = (sum(x)/N)
stdev = sum(x**2)/N

# print(x,y)
print('<x^2> = ', stdev)
print('<x>^2 = ', mean**2)
print('uncertainty = ',np.sqrt(stdev - mean**2))

plt.plot(x, y, 'o', color='lime')
plt.plot(x, y, '--', color='black')
plt.xticks(np.arange(-6, 6, 1))
plt.yticks(np.arange(-6, 6, 1))
s =  f'stdev = {round(stdev,2)}, \nmean^2 = {round(mean**2,2)},\nuncertainty = {round(np.sqrt(stdev - mean**2),2)}'
plt.annotate(s,[-5,3])
plt.show()

