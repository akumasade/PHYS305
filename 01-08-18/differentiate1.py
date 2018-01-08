import math
import matplotlib.pyplot as plt

def f(x):
    return math.cos(x)

def fp(x):
    return -math.sin(x)

x = 1.0
dx = 0.1

d = []
er = []

#want plot of abs(error) vs dx
while (dx >=1e-12):
    deriv = (f(x+dx)-f(x))/(dx)
    error = abs(deriv-fp(x))
    #print x, deriv, fp(x), abs(deriv-fp(x))
    d.append(dx)
    er.append(error)
    dx = dx*0.1

plt.figure()
plt.plot(d, er, 'ro')
plt.loglog()
plt.show()
