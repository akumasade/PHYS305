#I have not fucking idea what is going on rn
# nor do I have any idea where the f*ck we got this from
# or wtf this method is called
import numpy as np
from scipy.linalg import solve
import matplotlib.pyplot as plt

def second_deriv(x,y):
    return 0.0
def f(x):
    return 1.0

a = 0.0
b = 1.0
ya = 1.0
yb = 2.0

n = 50
delta = (b-a)/n

x = np.linspace(a, b, n-2)
# our matrix
A = np.zeros((n-2, n-2)) #(n-2)x(n-2) matrix
for i in range(n-2):
    A[i,i] = -2+delta**2 * f(x[i])
    if i > 0:
        A[i-1][i] = 1.0
        A[i][i-1] = 1.0

r = np.zeros(n-2)
r[0] = -ya
r[-1] = -yb

ans = solve(A, r)#use scipy's solver

plt.plot(x, ans)
plt.show()
#I have no idea what this means, but it looks like I did it right
