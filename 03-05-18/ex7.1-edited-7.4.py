import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.
# infinite square well with boundaries at -1 and +1

def deriv2(x, y):
    #return -y[0]
    return -y[2]*y[0] #<-- changed this

# Boundary conditions.
BC_a = -1.0
BC_ya = 0.0
BC_b = 1.0
BC_yb = 0.0

# Integration step and bisection tolerance.
dx = 0.01
tol = 1.e-4

# Range to search for a solution.
# intial guesses
zl = 0.0
zr = 5.0

#-----------------------------------------------------------------------

def f(x, y):
    return np.array([y[1], deriv2(x, y), 0.0])#<-- changed this

def rk4_step(func, x, y, dx):
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+0.5*dx, y+0.5*dy2)
    dy4 = dx*func(x+dx, y+dy3)
    y += (dy1 + 2*dy2 + 2*dy3 + dy4)/6.
    x += dx
    return x, y

def integrate(func, a, ya, z, b, dx):
    x = a
    y = np.array([ya, 1.0, z])#<---changed this
    while x < b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
    return x, y

#bisection fucntions
def g(z):
    x, y = integrate(f, BC_a, BC_ya, z, BC_b, dx)
    return y[0]-BC_yb

def bisect(func, zl, zr, tol):
    n = 0
    while zr-zl > tol:
	zm = 0.5*(zl + zr)
	if func(zm)*func(zl) > 0:
	    zl = zm
	else:
	    zr = zm
	n += 1
        print n, zl, func(zl), zr, func(zr)#func = g

    return n,zl,zr

#one secant iteration
def secant1(g, zl, zr):
    gl = g(zl)
    gr = g(zr)
    return zl + (zr-zl)*(-gl)/(gr-gl)

def main():

    global zl, zr
    n,zl,zr = bisect(g, zl, zr, tol)

    print "Root lies in range (%f, %f) after %d iterations"%(zl, zr, n)
    print "Function value =", g(0.5*(zl+zr))

    zs = secant1(g, zl, zr)
    print 'After secant iteration, z =', zs, 'function value =', g(zs)

    #plot the solution
    x = np.linspace(BC_a, BC_b)
    y = np.array([integrate(f, BC_a, BC_ya, zs, i, dx)[1][0] for i in x])
    plt.plot(x,y)
    plt.show()

if __name__ in ('__main__'):
    main()
