import numpy as np

#-------------part a--------------
def RungeKutta(f, xn, yn, dx):

    s = 6
    a = np.array([0., 1./5, 3./10, 3./5, 1., 7./8])
    b = np.array(
        [[0. for i in range(6)],
        [1./5, 0., 0., 0., 0., 0.],
        [3./40, 9./40, 0., 0., 0., 0.],
        [3./10, -9./10, 6./5, 0., 0., 0.],
        [-11./54, 5./2, -70./27, 35./27, 0., 0.],
        [1631./55296, 175./512, 575./13824, 44275./110592, 253./4096, 0.]])
    c = np.array([37./378, 0., 250./621, 125./594, 0., 512./1771])

    dy = np.zeros(s)
    for i in range(s):
        x = xn + a[i]*dx
        y = yn + dy.dot(b[i])
        print f(x, y)
        dy[i] = dx*f(x, y)

    return xn+dx, yn+dy.dot(c)

#from duffing.py
alpha = -2.0
beta = 1.0
delta = 0.0


iname = 'Runge Kutta'
xmax = 30.0
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

    x, y = RungeKutta(f, x, y, dx)

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
