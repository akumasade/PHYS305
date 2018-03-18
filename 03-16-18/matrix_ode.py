# Solve a simple boundary-value problem by matrix methods.

import math
import numpy as np
from scipy.linalg import solve, det
import matplotlib.pyplot as plt

def f(x):
    return 1.0

a = 0.0
b = 1.0
ya = 1.0
yb = 2.0

n = 100
delta = (b-a)/n
x = np.linspace(a, b, n)

# Simple difference matrix A (dimension n-2).

xx = x[1:-1]		# interior points
A = np.zeros((n-2,n-2))
for i in range(n-2):
    A[i,i] = -2.0 + delta**2*f(xx[i])
    if i > 0:
        A[i,i-1] = 1.0
    if i < n-3:
        A[i,i+1] = 1.0

print ''
print A
print np.linalg.det(A)
# Right-hand side r.
        
r = np.zeros(n-2)
r[0] = -ya
r[-1] = -yb

print r

y0 = solve(A, r)

# Compare the results with the analytic solution.

yy0 = np.array([ya])
yy0 = np.append(yy0, y0)
yy0 = np.append(yy0, yb)
A = (2-math.cos(1.0))/math.sin(1.0)
B = 1.0
yy1 = A*np.sin(x) + B*np.cos(x)

plt.plot(x, yy0, label='numerical')
plt.plot(x, yy1, label='analytic')
print 'maximum error =', np.max(np.abs(yy0-yy1))
plt.legend(loc='best')
plt.show()
