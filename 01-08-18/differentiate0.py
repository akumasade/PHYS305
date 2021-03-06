import math
import matplotlib.pyplot as plt

def f(x):
    return math.cos(x)

def fp(x):
    return -math.sin(x)

x = 1.0
dx = 1e-6
deriv = (f(x+dx)-f(x))/(dx)
print x, deriv, fp(x), abs(deriv-fp(x))
     
