import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.
# add in arguments for U0
a = 25.0
U0 = 0.0 #escape at energy = 0 #1 = 0 - (-1)

def U(x):
    return -np.exp(-np.sqrt(abs(x)))

def deriv2(x, y):
    return (-y[2] + U(x))*y[0]		# y[2] is the eigenvalue, z

# Integration step and bisection tolerance.

dx = 0.01
tol = 1.e-4


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
    #eta = math.sqrt(max(U0-z, 0.0))
    x, y = integrate(f, 0.0, y0, yp0, z, a, dx)
    #return y[1] + eta*y[0]
    return y[1] #goes to zero at x = a

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
    x = -a
    eta = 0.0 #math.sqrt(U0-z)
    y = np.array([1.0, eta, z])
    #y = np.array([y0, yp0, z])
    xx = [x]
    yy = [y[0]]
    while x < a-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
        xx.append(x)
        yy.append(y[0])
    if y0 == 1.0 :
        parity = " (even)"
    elif y0 == 0.0:
        parity = " (odd)"

    xx = np.array(xx)
    yy = np.array(yy)
    yy = yy/np.max(yy)#normalized
    '''
    #only plotting data between -5 and 5
    ind = np.where(abs(xx) <= 5.0)
    xx = xx[ind]
    yy = yy[ind]
    '''


    lab = 'z = {:.3f}'.format(z) + parity #label
    plt.plot(xx, yy, label= lab)
    #plt.xlim(-a, a)
    plt.xlabel('x')
    plt.ylabel('y')

def get_z():
    # Range to search for a solution.
    zlin = -1.0
    #zrin = U0 - 0.75

    zss = [] #array to hold all possible solutions (even)

    #loop thru all possible z
    zl = zlin
    zr = zl
    dz = 0.01
    while ( zl < 0.0):#only finding 5 sloutions for even/odd
        zr += dz
        if g(zr)*g(zl) <= 0.:
            n,zzl,zzr = bisect(g, zl, zr, 1.e-6)
            zs = secant1(g, zzl, zzr)
            print 'z =', zs, 'g(z) =', g(zs)
            zss.append(zs)
            plotz(f, zs)
        zl = zr
    return zss

def main():


    global y0, yp0, U0

    zsplot = []

    # find even solutions (y0 = 1, y'0 = 0 at x=0)
    y0 = 1.0
    yp0 = 0.0

    zss_even = get_z()


    # find odd solutions (y0 = 0, y'0 = 1 at x=0)
    y0 = 0.0
    yp0 = 1.0

    zss_odd = get_z()


    plt.title("Solutions")
    #plt.xlabel("x")
    #plt.ylabel("y")
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig("4-long.png", bbox_inches='tight')
    plt.show()



if __name__ in ('__main__'):
    main()
