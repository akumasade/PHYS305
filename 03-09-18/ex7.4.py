import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.

def deriv2(x, y):
    return -y[2]*y[0]				# y[2] is the eigenvalue z

# Boundary conditions.

BC_a = -1.0
BC_ya = 0.0
BC_b = 1.0
BC_yb = 0.0

# Integration step and bisection tolerance.

dx = 0.01
tol = 1.e-4

# Range to search for a solution.

zlin = 10.0
zrin = 30.0

#-----------------------------------------------------------------------

def f(x, y):
    return np.array([y[1], deriv2(x, y), 0.]) 	# z' = 0

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
    y = np.array([ya, 1.0, z])			# z is the eigenvalue
    while x < b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
    return x, y

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

    return n,zl,zr

def secant1(g, zl, zr):
    gl = g(zl)
    gr = g(zr)
    return zl + (zr-zl)*(-gl)/(gr-gl)

def plotz(f, z):
    x = BC_a
    y = np.array([BC_ya, 1.0, z])
    xx = [x]
    yy = [y[0]]
    while x < BC_b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
        xx.append(x)
        yy.append(y[0])
    plt.plot(xx, yy, label='z = {:.3f}'.format(z))
    plt.xlim(BC_a, BC_b)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.plot([BC_a], [BC_ya], 'o', clip_on=False,
             markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
    plt.plot([BC_b], [BC_yb], 'o', clip_on=False,
             markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')

def main():

    zl = zlin
    zr = zrin

    n,zzl,zzr = bisect(g, zl, zr, 1.e-6)
    #print "Root lies in (%f, %f) after %d iterations"%(zzl, zzr, n)
    #print "Function value =", g(0.5*(zzl+zzr))
    zs = secant1(g, zzl, zzr)
    print 'z =', zs, 'g(z) =', g(zs)
    plotz(f, zs)
    plt.legend(loc='best')

    plt.show()

if __name__ in ('__main__'):
    main()
