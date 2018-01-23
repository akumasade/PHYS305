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

xmax = 10.*math.pi
dx = 0.02

xplots = []
yplots = []
vplots = []
eplots = []

x = 0.0
y0 = 1.3
v0 = 1.5
y = np.array([y0, v0])
xplot = [x]
yplot = [y[0]]
vplot = [y[1]]
eplot = [energy(x, y)]
while x < xmax:
    yp = y[0]
    vp = y[1]
    if x == 0: yp = y0 + 1
    x, y = euler_step(f, x, y, dx)
    xplot.append(x)
    yplot.append(y[0])
    vplot.append(y[1])
    eplot.append(energy(x, y))
    if yp <= y0 and y[0] > y0:
        break

vint = vp + (y[1]-vp)*(y[0]-y0)/(y[0]-yp)
print 'vint =', vint, 'err =', abs((vint-v0)/v0)

plt.figure()
plt.plot(yplot, vplot)
plt.xlabel('y')
plt.ylabel('vy')
plt.show()
