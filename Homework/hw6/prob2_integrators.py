#The different integration functions used
#for problem 2
import numpy as np

def rk4_step(func, x, y, dx):
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+0.5*dx, y+0.5*dy2)
    dy4 = dx*func(x+dx, y+dy3)
    y += (dy1 + 2*dy2 + 2*dy3 + dy4)/6.
    x += dx
    return x, y

def midpoint_step(func, x, y, dx):
    dy = dx*func(x, y)
    y += dx*func(x+0.5*dx, y+0.5*dy)
    x += dx
    return x, y

def rk6_step(f, xn, yn, dx):

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

    dy = np.zeros((s,2))
    for i in range(s):
        x = xn + a[i]*dx
        y = yn + (dy*b[i].reshape(s,1)).sum(axis=0)
        dy[i] = dx*f(x, y)

    return xn+dx, yn+(dy*c.reshape(s,1)).sum(axis=0)
