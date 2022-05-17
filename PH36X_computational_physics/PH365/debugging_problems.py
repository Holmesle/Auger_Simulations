import numpy as np
import matplotlib.pyplot as plt

'''
Code 1 - range issue linspace only had one step, also missing a final parenthesis
'''

# Plot the function y(x) of slope 5 and y-intercept 3 from x = 0 to x = 10

# m = 5
# b = 3

# x_values = np.linspace(0, 10, 100)


# def y(x):
#     return m*x + b


# plt.plot(x_values, y(x_values))

# print('I am about to show the plot...')
# plt.show()

'''
Code 2 - missing a colon, inclusive requires n+1 because python stops before n ## double check, 
total = 0 needs to be outside the foor loop
'''
# print the sum of square roots from 1 to 4 inclusive

total=0

for x in range(1, 5):
    print(total, '+', np.sqrt(x))
    total= total + np.sqrt(x)

print(total)

'''
Code 3 - run time error, divide by zero
'''
# find the sum of values of 1/x from [-10, 10] (inclusive), excluding any undefined values

# total=0

# for x in range(-10, 10):
#     if x == 0:
#         total = total + 0
#     else:
#         total=total + 1/x

# print(total)

# '''
# Code 4 - I see nothing wrong with this
# '''

# plot the line y = 10x + 3

# def f(x):
#     m = 10
#     b = 3
#     return m*x + b

# x = np.linspace(1, 10, 10)
# plt.plot(f(x), x)
# plt.show()

'''
Code 5 - no asteric to signify multiplication, no plt.show(), arrange set up like linspace
'''

# Plot a parabola with zeros at 1 and -1 and a second derivative of 10

# y = 5(x-1)(x+1)

# x=np.arange(-5, 5, 0.1)

# def y(x):
#     a = (x - 1)
#     b = (x+1)
#     print(a, b)
#     c = a*b
#     return 5*c

# plt.plot(x, y(x))
# plt.show()

'''
Code 6 -  i not defined inside function, no change to i such that it reaches upper bound and exits
'''

# def factorial(n):
#     i = 1
#     product = 1
#     while i <= n:
#         product=product * i
#         i += 1
#     return product

# print('1! =', factorial(1))
# print('2! =', factorial(2))
# print('3! =', factorial(3))

'''
Code 7 - powers uses **, extra parenthesis
'''

# def hypotenuse(leg1, leg2):
#     return np.sqrt(leg1 ** 2 + leg2 ** 2)

# print('a triangle with legs 3 and 4 has hypotenuse', hypotenuse(3, 4))
# print('a triangle with legs pi and pi has hypotenuse',hypotenuse(np.pi, np.pi))
