import sys
import math
import numpy as np
import matplotlib.pyplot as plt
#-------------part a--------------
p = 2**(1./3)
q = 1.- p
r = 4 - 2*p


def kickdrift(f, x, y, dx):
    s = 4
    c = [1/r, q/r, q/r, 1/r]
    d = [0, 2/r, -2*p/r, 2/r]

    for i in range(s):
        y[1] = y[1] + d[i]*deriv2(x,y)*dx
        y[0] = y[0] + c[i]*y[1]*dx

    return x+dx, y

#-------------part b--------------
#from duffing.py
alpha = -2.0
beta = 1.0
delta = 0.0


iname = 'kick drift'
xmax = 20.0
y0 = 1.3
dx = 0.01

def deriv2(x, y):
    return -delta*y[1] - alpha*y[0] - beta*y[0]**3

def phi(x, y):
    return 0.5*alpha*y[0]**2 + 0.25*beta*y[0]**4

def energy(x, y):
    return 0.5*y[1]**2 + phi(x, y)

def f(x, y):
    return np.array([y[1], deriv2(x, y)])

print 'alpha =', alpha, 'beta =', beta
print 'integrator =', iname, 'dx =', dx

xplots = []
yplots = []
vplots = []
eplots = []

x = 0.0
y0 = 1.0
v0 = 1.5
y = np.array([y0, v0])
xplot = [x]
yplot = [y[0]]
vplot = [y[1]]
E0 = energy(x, y)
Escale = max(abs(E0), 1.0)
print 'E0, Escale =', E0, Escale
eplot = [0.0]
while x < xmax:

    x, y = kickdrift(f, x, y, dx)

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

plt.xlabel('x')
plt.ylabel('v')
#plt.xlim(-3.0, 3.0)
#plt.ylim(-2.5, 2.5)
plt.title('Duffing oscillator, integrator = '+iname)
subplot2 = fig.add_subplot(2, 1, 2)
for xplot,eplot in zip(xplots, eplots):
    plt.plot(xplot, eplot)
plt.xlabel('t')
plt.ylabel('$\delta$E/E$_0$')
plt.tight_layout()
plt.savefig("2b.png")
#c zxplt.show()

#-------------part b--------------
i= 0
dxs = [2**(-x) for x in range(1,14)]

xplots = []
yplots = []
vplots = []
eplots = []

diffs = []

for dx in dxs:
    x = 0.0
    y0 = 1.0
    v0 = 1.5
    y = np.array([y0, v0])
    xplot = [x]
    yplot = [y[0]]
    vplot = [y[1]]
    E0 = energy(x, y)
    Escale = max(abs(E0), 1.0)
    #print 'E0, Escale =', E0, Escale
    eplot = [0.0]
    while x < xmax:

        x, y = kickdrift(f, x, y, dx)

        xplot.append(x)
        yplot.append(y[0])
        vplot.append(y[1])
        eplot.append((energy(x, y)-E0)/Escale)
    E20 = energy(x,y)

    diffs.append(abs(E20-E0))

    xplots.append(xplot)
    yplots.append(yplot)
    vplots.append(vplot)
    eplots.append(eplot)

fig = plt.figure()

plt.plot(dxs, diffs, 'r')
plt.loglog()
plt.xlabel('$\Delta$t')
plt.ylabel('|E(20)-E(0)|')
plt.title('Duffing oscillator, integrator = '+iname)

plt.savefig("2c.png")
#plt.show()
#-------------part c--------------
i=0
xplots = []
yplots = []
vplots = []
eplots = []

#forward
x = 0.0
y0 = 1.0
v0 = 1.5
y = np.array([y0, v0])
xplot = [x]
yplot = [y[0]]
vplot = [y[1]]
E0 = energy(x, y)
Escale = max(abs(E0), 1.0)
#print 'E0, Escale =', E0, Escale
print "Initial Conditions:", x, y
eplot = [0.0]
while x < xmax:

    x, y = kickdrift(f, x, y, dx)

    xplot.append(x)
    yplot.append(y[0])
    vplot.append(y[1])
    eplot.append((energy(x, y)-E0)/Escale)

xplots.append(xplot)
yplots.append(yplot)
vplots.append(vplot)
eplots.append(eplot)

#reverse
dx = -dx


xplot = [x]
yplot = [y[0]]
vplot = [y[1]]
E0 = energy(x, y)
Escale = max(abs(E0), 1.0)
#print 'E0, Escale =', E0, Escale
eplot = [0.0]
while x >= 0.0:

    if x == 0.0: print "Reverse ending point:", x, y

    x, y = kickdrift(f, x, y, dx)

    xplot.append(x)
    yplot.append(y[0])
    vplot.append(y[1])
    eplot.append((energy(x, y)-E0)/Escale)

xplots.append(xplot)
yplots.append(yplot)
vplots.append(vplot)
eplots.append(eplot)


#It's time reversible~!! (though it probably should not be, goddammit)
