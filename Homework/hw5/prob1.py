import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.

def deriv2(x, y):
    #return -y[0]
    return -y[0] - 10*y[0]**3

# Boundary conditions.

BC_a = 0.0
BC_ya = 1.0
BC_b = 1.0
BC_yb = 2.0

# Integration step and bisection tolerance.

dx = 0.01
tol = 1.e-4

# Range to search for a solution.
# intial guesses
zl = 0.0
zr = 5.0
#gotta find all the roots now
#-----------------------------------------------------------------------

def f(x, y):
    return np.array([y[1], deriv2(x, y)])

def rk4_step(func, x, y, dx):
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+0.5*dx, y+0.5*dy2)
    dy4 = dx*func(x+dx, y+dy3)
    y += (dy1 + 2*dy2 + 2*dy3 + dy4)/6.
    x += dx
    return x, y

def integrate(func, a, ya, za, b, dx):
    x = a
    y = np.array([ya, za])
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

    #first need to plot g(z) to find the values

    x = np.linspace(-100, 100, 400)
    y = np.array([g(i) for i in x])
    zero = x*0.0
    plt.plot(x,y)
    plt.plot(x, zero)
    plt.xlabel("z")
    plt.ylabel("g(z)")
    plt.savefig("1-g_plot.png")
    plt.clf()

    #list of zl, zr to search thru (from plotting)
    zls = [-83.0, -48.0, -33.0, 13.0, 19.0, 49.0, 69.0]
    zrs = [-81.0, -42.0, -30.0, 16.0, 22.0, 53.0, 72.0]

    for l, r in zip(zls, zrs):
        n,l,r = bisect(g, l, r, tol)
        print "Root lies in range (%f, %f) after %d iterations"%(l, r, n)
        print "Function value =", g(0.5*(l+r))
        zs = secant1(g, l, r)
        print 'After secant iteration, z =', zs, 'function value =', g(zs)
        #plotting
        x = np.linspace(BC_a, BC_b, 100)
        y = np.array([integrate(f, BC_a, BC_ya, zs, i, dx)[1][0] for i in x])
        plt.plot(x,y, label="z = %s"%zs)
        plt.plot(x[0], y[0], 'ro')
        plt.plot(x[-1], y[-1], 'ro')

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig("1-solutions.png", bbox_inches='tight')
    plt.show()
    
if __name__ in ('__main__'):
    main()
