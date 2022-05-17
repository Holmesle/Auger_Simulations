import numpy as np
import matplotlib.pylab as plt


# 09/10/2021 - 06:01PM to 7:35PM - Leah Holmes

'''
This program construct a triangle out of three squares made up of units, 
such that the perimeter of the triangle is 1000 units.

'''
# requested perimeter 

perimeter = 1000

'''
1. Set up original squares (5x4x3 center triangle)
'''

# define sides of blocks [BL->TL, TL->TR, TR->BR, BR->BL]
b1 = np.array(([[1, 1], [1, 6]], [[1, 6], [6, 6]],
              [[6, 6], [6, 1]], [[6, 1], [1, 1]]))
b2 = np.array(([[1, 1], [1, 5]], [[1, 5], [5, 5]],
              [[5, 5], [5, 1]], [[5, 1], [1, 1]]))
b3 = np.array(([[1, 1], [1, 4]], [[1, 4], [4, 4]],
              [[4, 4], [4, 1]], [[4, 1], [1, 1]]))

'''
2. Create transformation funcitons
'''
# function for rotating blocks

def Rotate(block,deg):
    theta = deg*np.pi/180
    final = []
    rot = np.array([[np.cos(theta), -np.sin(theta)],
                    [np.sin(theta), np.cos(theta)]])
    for v in block: 
        final.append(np.dot(rot,v))
    return np.array(final)

# function for translating blocks

def Translate(block,x,y,scale):
    final = []
    shift = np.array([[x, x], [y, y]])
    for v in block: 
        final.append(v + shift)
    return np.array(final)

def Scale(block,h,w):
    final = []
    k = np.array([[h,0],[0,w]])
    for v in block:
        final.append(np.dot(k, v))
    return np.array(final)

# apply transforms to obtain original

b1 = Translate(b1, 0, 3.3, 0) # blue
b2 = Rotate(b2, 54)  # red
b2 = Translate(b2, 9.48, 0.55, 0) # red
b3 = Rotate(b3, 54)  # lime
b3 = Translate(b3, 8.65, 6.15, 0) # lime

# scale original by 1000/12 to obtain a perimeter of 1000

k = perimeter/12
b1_scaled = Scale(b1, k, k)
b2_scaled = Scale(b2, k, k)
b3_scaled = Scale(b3, k, k)

'''
Check length of perimeter
'''
def Perimeter(block1,block2,block3):
    side1_length = np.sqrt((block1[1][0][1] - block1[1][0][0])**2 
                        + (block1[1][1][1] - block1[1][1][0])**2)
    side2_length = np.sqrt((block2[2][0][1] - block2[2][0][0])**2
                       + (block2[2][1][1] - block2[2][1][0])**2)
    side3_length = np.sqrt((block3[0][0][1] - block3[0][0][0])**2
                       + (block3[0][1][1] - block3[0][1][0])**2)
    return round(side1_length + side2_length + side3_length,3)

length = Perimeter(b1, b2, b3)
length_scaled = Perimeter(b1_scaled, b2_scaled, b3_scaled)

print('The perimeter of the center triangle is %.2f' %length_scaled)
'''
4. plot squares to check orientation and display
'''

# function for graphing a single given block

def graphblock(block, clr, ax):
    for v in block:
        ax.plot(v[0], v[1], color=clr)


# # graph blocks original
# fig, ax1 = plt.subplots()
# graphblock(b1, 'blue', ax1)
# graphblock(b2, 'red', ax1)
# graphblock(b3, 'lime', ax1)
# ax1.set_aspect('equal', adjustable = 'box')
# ax1.set_title(f'3 Blocks (perimeter = %{length})')

# # graph blocks scaled
fig1, ax2 = plt.subplots()
graphblock(b1_scaled, 'blue', ax2)
graphblock(b2_scaled, 'red', ax2)
graphblock(b3_scaled, 'lime', ax2)
ax2.set_aspect('equal', adjustable='box')
ax2.set_title(f'3 Blocks (perimeter = {length_scaled})')
plt.show()
