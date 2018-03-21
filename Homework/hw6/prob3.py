# Based off of solution to hw5.4

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.

def deriv2(x, y):
    if x==0:
        return 0.0 #division by zero
    else:
        return -(y[2] + 2.0/x)*y[0]		# y[2] is the eigenvalue z

# Integration step and bisection tolerance.

dx = 0.005
tol = 1.e-8

# Range to search for a solution.

zlin = -10.0
zrin = -1.e-8

x0 = 30.0		# integration to too large x0 becomes unstable

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

def integrate(func, a, ya, ypa, z, b, dx):
    x = a
    y = np.array([ya, ypa, z])			# z is the eigenvalue
    while x < b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
    return x, y

def g(z):
    x, y = integrate(f, 0.0, 1.-parity, parity, z, x0, dx)
    eta = math.sqrt(-z)
    return y[1]+eta*y[0]			# match to exterior solution

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
    x = 0.0
    y = np.array([1.-parity, parity, z])
    xx = [x]
    yy = [y[0]]
    while x < x0-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
        xx.append(x)
        yy.append(y[0])
    plt.plot(xx, yy, label='z = {:.3f}'.format(z))
    plt.xlim(0., x0)
    plt.xlabel('x')
    plt.ylabel('y')

def main():
    global parity

    #plotting g(z)
    plt.figure(figsize=(12,5))

    ng = 100
    zz = 10.**np.linspace(math.log10(-zlin), math.log10(-zrin), ng)
    gge = np.zeros(ng)
    ggo = np.zeros(ng)
    for i in range(ng):
        parity = 0
        gge[i] = g(-zz[i])
        parity = 1
        ggo[i] = g(-zz[i])

    plt.subplot(1,2,1)
    plt.semilogx(zz, gge, label='even')
    plt.semilogx(zz, ggo, label='odd')
    plt.semilogx([-zlin,-zrin], [0.,0.], 'r--')
    plt.xlim(-zlin, -zrin)
    plt.xlabel('-z')
    plt.ylabel('g(z)')
    plt.ylim(-1.0, 1.0)
    plt.legend(loc='best')

    # Use gge and ggo to bracket the roots.

    plt.subplot(1,2,2)
    
    print "Even"
    parity = 0
    for i in range(ng-1):
        if gge[i]*gge[i+1] <= 0:
            nn,zzl,zzr = bisect(g, -zz[i], -zz[i+1], 1.e-6)
            zs = secant1(g, zzl, zzr)
            print 'z =', zs, 'g(z) =', g(zs)
            if zs <= -0.10: #plotting only first 5 solutions
                plotz(f, zs)
    
    print "Odd"
    parity = 1
    for i in range(ng-1):
        if ggo[i]*ggo[i+1] <= 0:
            nn,zzl,zzr = bisect(g, -zz[i], -zz[i+1], 1.e-6)
            zs = secant1(g, zzl, zzr)
            print 'z =', zs, 'g(z) =', g(zs)
            if zs <= -0.10: #plotting only first 5 solutions
                plotz(f, zs)

    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig('3b.png', bbox_inches='tight')
    plt.show()
    
    
    
if __name__ in ('__main__'):
    main()