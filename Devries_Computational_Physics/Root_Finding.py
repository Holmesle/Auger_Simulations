from xmlrpc.server import DocXMLRPCRequestHandler
import numpy as np
import matplotlib.pyplot as plt

def bisection(f,xi):
    # find the root of the equation using bisection method
    a = f[0]
    b = f[-1]
    if f(xi) < 0:
        pass

print(np.pi)