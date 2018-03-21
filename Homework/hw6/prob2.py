#Based off of ex7.4 solution
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from prob2_integrators import *


#-----------------------------------------------------------------------
#
# Define the problem.
def deriv2(x, y):
    return -(2.0*y[1] + y[0] - 10.0*y[0]**2. + 5.0*y[0]**3.)# y[2] is the eigenvalue, z

# Boundary conditions.

BC_a = 0.0
BC_ya = 0.0
BC_b = 2.0
BC_yb = 1.0

# Integration step and bisection tolerance.
#dx = 0.01 #moved to main
tol = 1.e-4

# Range to search for a solution.
# assumed to have solution between 60 and 70
zlin = 60.0
zrin = 70.0

#-----------------------------------------------------------------------

def f(x, y):
    return np.array([y[1], deriv2(x, y)]) 	# z' = 0


def integrate(func, a, ya, ypa, b, dx):

    x = a
    y = np.array([ya, ypa])			# z is the eigenvalue
    while x < b-0.5*dx:
        x, y = integrator(f, x, y, dx)#perform a time step of size dx
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
    y = np.array([BC_ya, z])
    xx = [x]
    yy = [y[0]]
    while x < BC_b-0.5*dx:
        x, y = integrator(f, x, y, dx)
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
    
def get_zs():
    zl = zlin
    zr = zrin

    n,zzl,zzr = bisect(g, zl, zr, 1.e-6)
    zs = secant1(g, zzl, zzr)
    print 'z =', zs, 'g(z) =', g(zs)
    #plotz(f, zs)
    return zs

def main():

    global integrator, dx
    
    #variables to keep track once we setup the loop
    integrators = [midpoint_step, rk4_step, rk6_step]
    dxs = [0.1, 0.01, 0.001]
    max_errs = np.zeros((3,3))#array to hold our table of max errors
    zss = np.zeros((3,3))#array to hold our solutions
    
    
    '''
    #plotting g(z)
    x = np.linspace(60, 70, 100)
    gz = [g(i) for i in x]


    plt.plot(x, gz)
    plt.plot(x, 0*x)
    plt.xlabel('z')
    plt.ylabel('g(z)')
    #plt.ylim(-1.0, 1.0)
    plt.legend(loc='best')
    plt.show()
    '''
    #get best solution
    integrator = rk6_step
    dx = 0.001
    zbest = get_zs()
    #integrate
    x = BC_a
    y = np.array([BC_ya, zbest])
    xx = [x]
    ybest = [y[0]]#best solution
    while x < BC_b-0.5*dx:
        x, y = integrator(f, x, y, dx)
        xx.append(x)
        ybest.append(y[0])
    
    
    
    
    
    #find the difference between them and the "most accurate solution"
    for i,igt in enumerate(integrators):
        for j,d in enumerate(dxs):
            dx = d
            integrator = igt
            z = get_zs() #find the solution
            zss[i][j] = z
            
            if j==0: n = 100#dx = 0.1   
            elif j==1: n = 10#dx = 0.01   
            elif j==2: n = 1#dx=0.001
            yi = ybest[0::n]
            
            #perform plotting (y - y_best)
            x = BC_a
            y = np.array([BC_ya, z])
            xx = [x]
            ind  = 0
            yy = [y[0] - yi[ind]]#current solution - best solution
            
            while x < BC_b-0.5*dx:
                x, y = integrator(f, x, y, dx)
                ind += 1
                xx.append(x)
                yy.append(y[0]- yi[ind])#find the difference
            
            #calculate maximum absolute difference
            yarray = np.array(yy)
            max_errs[i][j] = np.max(abs(yarray))
            
            #plotting
            if i==0: int_name = "Midpoint"  
            elif i==1: int_name = "Runge-Kutta 4"  
            elif i==2: int_name = "Runge-Kutta 6"
            
            plt.plot(xx, yy, label='%s, dx = %s, z = %s'%(int_name, dx, z))
            plt.xlim(BC_a, BC_b)
            plt.xlabel('x')
            plt.ylabel('y')
            
            #plt.plot([BC_a], [BC_ya], 'o', clip_on=False,
            #         markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
            #plt.plot([BC_b], [BC_yb], 'o', clip_on=False,
            #         markeredgewidth=1, markeredgecolor='r', markerfacecolor='None')
            
            
    print "Maximum Absolute Difference"
    print max_errs

    plt.title("Deviation from Runge-Kutta 6 , dx = 0.001")
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig("2.png", bbox_inches='tight')
    plt.show()

if __name__ in ('__main__'):
    main()
