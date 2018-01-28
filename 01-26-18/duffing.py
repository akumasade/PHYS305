import sys
import math
import numpy as np
import matplotlib.pyplot as plt

alpha = -2.0
beta = 1.0
delta = 0.0

def deriv2(x, y):
    return -delta*y[1] - alpha*y[0] - beta*y[0]**3

def phi(x, y):
    return 0.5*alpha*y[0]**2 + 0.25*beta*y[0]**4

def energy(x, y):
    return 0.5*y[1]**2 + phi(x, y)

def f(x, y):
    return np.array([y[1], deriv2(x, y)])

#----------------------------------------------------------------------

# Integrators (names are self-descriptive):

def implicit_euler_step(func, x, y, dx):	# 0
    dy = dx*func(x, y)
    y[1] += dy[1]
    y[0] += dx*y[1]
    x += dx
    return x, y

def euler_step(func, x, y, dx):			# 1
    y += dx*func(x, y)
    x += dx
    return x, y

def midpoint_step(func, x, y, dx):		# 2
    dy = dx*func(x, y)
    y += dx*func(x+0.5*dx, y+0.5*dy)
    x += dx
    return x, y

def rk3_step(func, x, y, dx):			# 3
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+dx, y-dy1+2*dy2)
    y += (dy1 + 4*dy2 + dy3)/6.
    x += dx
    return x, y

def rk4_step(func, x, y, dx):			# 4
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+0.5*dx, y+0.5*dy2)
    dy4 = dx*func(x+dx, y+dy3)
    y += (dy1 + 2*dy2 + 2*dy3 + dy4)/6.
    x += dx
    return x, y

#----------------------------------------------------------------------

def main(alpha_in, beta_in, dx, integrator, xmax, y0):

    alpha = alpha_in
    beta = beta_in

    if integrator == 0:
        iname = 'implicit Euler'
    elif integrator == 1:
        iname = 'Euler'
    elif integrator == 2:
        iname = 'Mid-point'
    elif integrator == 3:
        iname = 'Runge-Kutta 3'
    else:
        iname = 'Runge-Kutta 4'

    print 'alpha =', alpha, 'beta =', beta
    print 'integrator =', iname, 'dx =', dx

    xplots = []
    yplots = []
    vplots = []
    eplots = []

    for v0 in [0.25, 0.5, 1.0, 1.25, 1.5, 2.0, 3.0]:
        x = 0.0
        y = np.array([y0, v0])
        xplot = [x]
        yplot = [y[0]]
        vplot = [y[1]]
        E0 = energy(x, y)
        Escale = max(abs(E0), 1.0)
        print 'E0, Escale =', E0, Escale
        eplot = [0.0]
        while x < xmax:
            yp = y[0]
            vp = y[1]
            if x == 0: yp = y0 + 1
            if integrator == 0:
                x, y = implicit_euler_step(f, x, y, dx)
            elif integrator == 1:
                x, y = euler_step(f, x, y, dx)
            elif integrator == 2:
                x, y = midpoint_step(f, x, y, dx)
            elif integrator == 3:
                x, y = rk3_step(f, x, y, dx)
            else:
                x, y = rk4_step(f, x, y, dx)

            xplot.append(x)
            yplot.append(y[0])
            vplot.append(y[1])
            eplot.append((energy(x, y)-E0)/Escale)

        xplots.append(xplot)
        yplots.append(yplot)
        vplots.append(vplot)
        eplots.append(eplot)

    fig = plt.figure(figsize=(6,9))
    subplot1 = fig.add_subplot(2, 1, 1)
    for yplot,vplot in zip(yplots, vplots):
        plt.plot(yplot, vplot)
    if alpha < 0:
        yeq = math.sqrt(-alpha/beta)
        plt.plot(yeq, 0., 'kx')
        plt.plot(-yeq, 0., 'kx')
    plt.xlabel('y')
    plt.ylabel('v$_y$')
    plt.xlim(-3.0, 3.0)
    plt.ylim(-2.5, 2.5)
    plt.title('Duffing oscillator, integrator = '+iname)
    subplot2 = fig.add_subplot(2, 1, 2)
    for xplot,eplot in zip(xplots, eplots):
        plt.plot(xplot, eplot)
    plt.xlabel('x')
    plt.ylabel('$\delta$E/E$_0$')
    plt.tight_layout()
    plt.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
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
    result.add_option("-i",
                      dest="integrator", type="int", default ="2",
                      help="integrator [%default]")
    result.add_option("-x",
                      dest="xmax", type="float", default ="30.0",
                      help="integration interval [%default]")
    result.add_option("-y",
                      dest="y0", type="float", default ="1.3",
                      help="integration interval [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)#pass the args thru main function
