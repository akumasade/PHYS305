import math
import matplotlib.pyplot as plt

def f(x):
    return math.cos(x)

def fi(x):
    return math.sin(x)

#error should look like line
def trap(func, xmin, xmax, n):
    sum = 0.5*(f(xmin)+f(xmax))
    x = xmin
    dx = (xmax-xmin)/n
    for i in range(1,n):
        x = xmin + i*dx
        sum += f(x)#trapezoidal method
    return dx*sum

xmin = 0.0
xmax = 10.0
n = 4
nmax = 1000000
dxplot = []
errplot = []
while n < nmax:
    n *= 2
    #print n
    dx = (xmax-xmin)/n
    dxplot.append(dx)
    num_int = trap(f, xmin, xmax, n)
    true_int = fi(xmax)-fi(xmin)
    errplot.append(abs(num_int - true_int))

plt.figure()
plt.plot(dxplot, errplot)
plt.loglog()
plt.xlabel('$\Delta$x')
plt.ylabel('error')
plt.show()
