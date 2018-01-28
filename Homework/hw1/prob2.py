#Problem 2
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.cos(x)*x

def fi(x):
    return x*np.sin(x)+np.cos(x)

def basic(func, xmin, xmax, n):
    sum = 0.
    x = xmin
    dx = (xmax-xmin)/n
    for i in range(n+1):
        x = xmin + i*dx
        sum += f(x)
    return dx*sum

def trap(func, xmin, xmax, n):
    sum = 0.5*(f(xmin)+f(xmax))
    x = xmin
    dx = (xmax-xmin)/n
    for i in range(1,n):
        x = xmin + i*dx
        sum += f(x)#trapezoidal method
    return dx*sum

def simp(func, xmin, xmax, n):
    sum = f(xmin)+f(xmax)
    #sum = 0.
    x = xmin
    dx = (xmax-xmin)/n
    for i in range(1, n):
        x = xmin + i*dx
        if(i%2==0): sum += 2*f(x)
        else: sum += 4*f(x)
    return dx*sum/3.0

xmin = 0.0
xmax = 2.0
n = 4
nmax = 2**20
dxplot = []
errplot = []
er_basic = []
er_trap = []
er_simp = []

while n < nmax:
    true_int = fi(xmax)-fi(xmin)
    #dx = (xmax-xmin)/nn
    dx = 2.0/n
    dxplot.append(dx)

    er_basic.append(abs(basic(f, xmin, xmax, n)-true_int))
    er_trap.append(abs(trap(f, xmin, xmax, n)-true_int))
    er_simp.append(abs(simp(f, xmin, xmax, n)-true_int))

    n *= 2

plt.figure()
plt.plot(dxplot, er_basic, label="basic")
plt.plot(dxplot, er_trap, label="trapezoidal")
plt.plot(dxplot, er_simp, label="simpson")
plt.loglog()
plt.title("Problem 2a.")
plt.xlabel('$\Delta$x')
plt.ylabel('error')
plt.legend()
plt.savefig('prob2a.png')
#plt.show()
plt.clf()

#-------part b-------
Ns = np.array([4, 8, 16, 32])
dxs = Ns/2.0
traps = np.array([])
for n in Ns:
    traps = np.append(traps, [trap(f, xmin, xmax, n)])

print len(dxs), len(traps)
p = np.polyfit(dxs, traps, 3)
print p

dxx = np.linspace(-50, 50)
der = p[0]*dxx**3 + p[1]*dxx**2 + p[2]*dx + p[3]

actual = fi(xmax) - fi(xmin)
print "Error = ",abs(p[3]-actual)#print the dx=0 value

plt.plot(dxx, der, '-', label="answer at (0, %s)"%(p[3]))
plt.title('Problem 2b.')
plt.xlabel('$\Delta$x')
plt.ylabel('Integral')
plt.legend()
#plt.savefig('prob1c.png')
plt.show()
