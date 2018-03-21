# Based off of solutions to Homework 2, problem 2
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

alpha = -2.0
beta = 1.0
delta = 0.0

# Note: deriv2 is acc here, independent of x (= time).

def deriv2(y):
    return -delta*y[1] - alpha*y[0] - beta*y[0]**3

def phi(x, y):
    return 0.5*alpha*y[0]**2 + 0.25*beta*y[0]**4

def energy(x, y):
    return 0.5*y[1]**2 + phi(x, y)

def f(x, y):
    return np.array([y[1], deriv2(y)])

#----------------------------------------------------------------------

# integrator
w = np.array([-1.17767998417887, 0.235573213359357, 0.784513610477560])
wsum = 1.0 - 2.0*np.sum(w)

c = [0.5*w[2], 0.5*(w[2] + w[1]), 0.5*(w[1] + w[0]), 0.5*(w[0] + wsum)]
c = np.array(c + list(reversed(c))) #second half of c looks like fist half in reverse

d = np.array([0.0, w[2], w[1], w[0], wsum, w[0], w[1], w[2]])

def symp_step(func, x, y, dx):
    for cc,dd in zip(dx*c, dx*d):
        y[1] += dd*func(x, y)[1]
        y[0] += cc*y[1]
    x += dx
    return x, y

#----------------------------------------------------------------------

def main(alpha_in, beta_in, dx, xmax):

    alpha = alpha_in
    beta = beta_in
    print 'alpha =', alpha, 'beta =', beta

    x0 = 0.0
    y0 = 1.0
    yp0 = 1.5
    dx0 = dx

    # Part (c):

    print ''
    dxplot = []
    deplot = []
    #for  n in range(1,14):
    for  n in range(1,8):
        dx = 2.**(-n)
        print '(c) dx =', dx
        x = x0
        y = np.array([y0, yp0])
        E0 = energy(x, y)
        while x < xmax-0.5*dx:
            x, y = symp_step(f, x, y, dx)
        dxplot.append(dx)
        deplot.append(abs(energy(x, y)-E0))

    plt.loglog(dxplot, deplot, 'r+')
    plt.xlabel('$\delta$ x')
    plt.ylabel('$\delta$ E')
    #plt.title('Homework 2.2(c)')

    logdx = np.log10(np.array(dxplot[1:9]))
    logde = np.log10(np.array(deplot[1:9]))
    p = np.polyfit(logdx, logde, 1)
    print '(c) slope =', p[0]

    plt.tight_layout()
    plt.savefig("1a.png")
    plt.show()
    #plt.cla()

    # Part (d):

    x = x0
    y = np.array([y0, yp0])
    dx = dx0

    # Forward integration:

    while x < xmax-0.5*dx:
        x, y = symp_step(f, x, y, dx)

    # Backward integration:

    while x > 0.5*dx:
        x, y = symp_step(f, x, y, -dx)

    print '\n(d)', x, y
    print '(d)', y[0]-y0, y[1]-yp0
    tol = 5.e-13
    if math.sqrt(((y[0]-y0)**2).sum()) < tol \
         and math.sqrt(((y[1]-yp0)**2).sum()) < tol:
        print '(d) time reversible'
    else:
        print '(d) not time reversible'

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-a",
                      dest="alpha_in", type="float", default ="-2.0",
                      help="alpha_in [%default]")
    result.add_option("-b",
                      dest="beta_in", type="float", default ="1.0",
                      help="beta_in [%default]")
    result.add_option("-d",
                      dest="dx", type="float", default ="0.01",
                      help="step length [%default]")
    result.add_option("-x",
                      dest="xmax", type="float", default ="20.0",
                      help="integration interval [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
