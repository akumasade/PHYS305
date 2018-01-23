import math
import numpy as np
import matplotlib.pyplot as plt

delta = 0.0
alpha = -2.0
beta = 1.0

def deriv2(x, y):
    return -delta*y[1] - alpha*y[0] - beta*y[0]**3

def phi(x, y):
    return 0.5*alpha*y[0]**2 + 0.25*beta*y[0]**4

def energy(x, y):
    return 0.5*y[1]**2 + phi(x, y)

def f(x, y):
    return np.array([y[1], deriv2(x, y)])

def euler_step(func, x, y, dx):
    y += dx*func(x, y)
    x += dx
    return x, y

def implicit_euler_step(func, x, y, dx):
    dy = dx*func(x, y)
    x += dx
    y[1] += dy[1]
    y[0] += dx*func(x, y)[0]
    return x, y

def midpoint_step(func, x, y, dx):
    dy = dx*func(x, y)
    y += dx*func(x+0.5*dx, y+0.5*dy)
    x += dx
    return x, y

xmax = 10.*math.pi
dx = 0.025#we see that error is linear (proportinal to dx for explicit euler)

#loop over different initial velocities, x
initialxs = [float(x)*np.pi for x in range(11)]
initialvs = [0.25, 0.5, 1.0, 1.25, 1.5, 2.0]#4.2
initialvs = initialvs + [1.35, 1.355, 1.36]#4.3

#initialize plot
plt.figure()
for xinit in initialxs:
    for vinit in initialvs:

        xplots = []
        yplots = []
        vplots = []
        eplots = []

        x = xinit
        y0 = 1.3
        v0 = vinit
        y = np.array([y0, v0])
        xplot = [x]
        yplot = [y[0]]
        vplot = [y[1]]
        eplot = [energy(x, y)]
        while x < xmax:
            yp = y[0]
            vp = y[1]
            if x == 0: yp = y0 + 1
            #x, y = euler_step(f, x, y, dx)
            #x, y = implicit_euler_step(f, x, y, dx)
            x, y = midpoint_step(f, x, y, dx)
            xplot.append(x)
            yplot.append(y[0])
            vplot.append(y[1])
            eplot.append(energy(x, y))
            if yp <= y0 and y[0] > y0: break

        #vint = vp + (y[1]-vp)*(y[0]-y0)/(y[0]-yp)
        #print 'vint =', vint, 'err =', abs((vint-v0)/v0)

        plt.plot(yplot, vplot)#plot

plt.xlabel('y')
plt.ylabel('vy')
plt.show()
