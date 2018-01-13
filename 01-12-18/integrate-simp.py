import math
import matplotlib.pyplot as plt

def f(x):
    return math.cos(x)

def fi(x):
    return math.sin(x)

#error should look like line with tail on left
def simpson(func, xmin, xmax, n):
    sum = f(xmin)+f(xmax)
    #sum = 0.
    x = xmin
    dx = (xmax-xmin)/n
    for i in range(1, n):
        x = xmin + i*dx
        if(i%2==0): sum += 2*f(x)
        else: sum += 4*f(x)
    return dx*sum/3.0

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
    num_int = simpson(f, xmin, xmax, n)
    true_int = fi(xmax)-fi(xmin)
    errplot.append(abs(num_int - true_int))

plt.figure()
plt.plot(dxplot, errplot)
plt.loglog()
plt.xlabel('$\Delta$x')
plt.ylabel('error')
plt.show()
