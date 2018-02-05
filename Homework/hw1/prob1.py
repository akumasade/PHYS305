import math
import matplotlib.pyplot as plt
import numpy as np

def f(x):
    return np.exp(-x)*x*x#function: x^2 * e^(-x)

def fp(x):
    return -np.exp(-x)*x*(x-2)#actual derivative from wolfram alpha

def center(x, dx):
    return (f(x+dx)-f(x-dx))/(2*dx)

x = 3.0#point at which we're taking the derivative
dx = 0.5

d = []
er1 = []
er2 = []

#want plot of abs(error) vs dx
while (dx >=2**(-50)):
    onesided = (f(x+dx)-f(x))/(dx)#1
    centered = (f(x+dx)-f(x-dx))/(2*dx)#2
    error_1 = abs(onesided-fp(x))
    error_2 = abs(centered-fp(x))
    d.append(dx)
    er1.append(error_1)
    er2.append(error_2)
    dx /= 2

plt.figure()
plt.plot(d, er1, 'r-', label="onesided $e^{(1)}$")
plt.plot(d, er2, 'b-', label="centered $e^{(2)}$")
plt.loglog()
plt.title('Problem 1a.')
plt.xlabel('$\Delta$x')
plt.ylabel('Error')
plt.legend()
plt.savefig('prob1a.png')
#plt.show()
plt.clf()

#-------part b-------
dxs = np.array([0.1, 0.05, 0.025])
deriv  = center(x,dxs)
p = np.polyfit(dxs, deriv, 2)
#print p #[  2.49025873e-02  -5.81800652e-07  -1.49361196e-01]

dxx = np.linspace(-50, 50)
der = p[0]*dxx*dxx + p[1]*dxx + p[2]
#when dx=0, deriv=-1.49361196e-01
print p[2]
print "Error = ",abs(p[2]-fp(x))

plt.plot(dxx, der, 'g-', label="answer at (0, %s)"%(p[2]))
plt.title('Problem 1b.')
plt.xlabel('dx')
plt.ylabel('Derivative')
plt.legend()
plt.savefig('prob1b.png')
#plt.show()
plt.clf()
#-------part c-------
dxs = np.array([0.1, 0.05, 0.025, 0.0125])
deriv  = center(x,dxs)
p = np.polyfit(dxs, deriv, 3)
#print p

dxx = np.linspace(-50, 50)
der = p[0]*dxx**3 + p[1]*dxx**2 + p[2]*dx + p[3]

print p[3]
print "Error = ",abs(p[3]-fp(x))#print the dx=0 value

plt.plot(dxx, der, '-', label="answer at (0, %s)"%(p[3]))
plt.title('Problem 1c.')
plt.xlabel('$\Delta$x')
plt.ylabel('Derivative')
plt.legend()
plt.savefig('prob1c.png')
#plt.show()
