import sys
import math
import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------------------
#
# Define the problem.
# infinite square well with boundaries at -1 and +1

#x0 = 6.0
x0 = 4.0
def deriv2(x, y):
    #return -y[0]
    return (-y[2] + U(x))*y[0] #<-- changed this

#"confined" harmonic oscillator
def U(x):
    return x*x

# Integration step and bisection tolerance.
dx = 0.01
tol = 1.e-4


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

def integrate(func, a, ya, ypa, z, b, dx):
    x = a
    y = np.array([ya, ypa, z])			# z is the eigenvalue
    while x < b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
    return x, y

#bisection fucntions
def g(z):
    x, y = integrate(f, 0.0, y0, yp0, z, x0, dx)
    #print "y(x0) = ", y[0]
    return y[0]  #0 at x0
    #return y[1]#slope goes to zero at edge

def bisect(func, zl, zr, tol):
    n = 0
    while zr-zl > tol:
	zm = 0.5*(zl + zr)
	if func(zm)*func(zl) > 0:
	    zl = zm
	else:
	    zr = zm
	n += 1
        #print n, zl, func(zl), zr, func(zr)#func = g

    return n,zl,zr

#one secant iteration
def secant1(g, zl, zr):
    gl = g(zl)
    gr = g(zr)
    return zl + (zr-zl)*(-gl)/(gr-gl)

def plotz(f, z):
    a = x0
    x = -a
    #eta = math.sqrt(U0-z)
    eta = 0.0
    y = np.array([1.0, eta, z])
    #yl = np.exp(eta*x) #y value at left edge
    #y = np.array([yl, eta, z])
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

    lab = 'z = {:.3f}'.format(z) + parity
    plt.plot(xx, yy, label= lab)
    #plt.xlim(-a, a)
    plt.xlabel('x')
    plt.ylabel('y')

def get_z():
    zss = [] #array to hold all possible solutions (even)

    #loop thru all possible z
    zl = 0.0
    zr = zl
    dz = 0.25
    while len(zss) < 5:
        zr += dz
        #print zl, zr
        if g(zr)*g(zl) <= 0.:
            n,zzl,zzr = bisect(g, zl, zr, 1.e-6)
            zs = secant1(g, zzl, zzr)
            print 'z =', zs, 'g(z) =', g(zs)
            zss.append(zs)
            plotz(f, zs)
        zl = zr
    return zss


def main():

    global  y0, yp0


    # find even solutions (y0 = 1, y'0 = 0 at x=0)
    y0 = 1.0
    yp0 = 0.0

    zss_even = get_z()

    # find odd solutions (y0 = 0, y'0 = 1 at x=0)
    y0 = 0.0
    yp0 = 1.0

    zss_odd = get_z()


    #plt.xlabel("x")
    #plt.ylabel("y")
    plt.title("$x_{0} = $%s"%x0)
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig("3-solutions-%s.png"%x0, bbox_inches='tight')
    plt.show()


if __name__ in ('__main__'):
    main()
