i = 1


def factorial(n):
    product = 1
    for i in range(n+1):
        product = product * i
    print('inside i is', i)
    return product


print('i is', i)
print('3! is', factorial(3))
print('i is', i)
