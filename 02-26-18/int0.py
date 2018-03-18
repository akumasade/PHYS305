import numpy as np

def deriv2(x, y):			# differential equation to solve
    return -y[0]

def f(x, y):
    return np.array([y[1], deriv2(x, y)])

def rk4_step(func, x, y, dx):
    dy1 = dx*func(x, y)
    dy2 = dx*func(x+0.5*dx, y+0.5*dy1)
    dy3 = dx*func(x+0.5*dx, y+0.5*dy2)
    dy4 = dx*func(x+dx, y+dy3)
    y += (dy1 + 2*dy2 + 2*dy3 + dy4)/6.
    x += dx
    return x, y

def integrate(func, a, ya, za, b, dx):
    x = a
    y = np.array([ya, za])
    while x < b-0.5*dx:
        x, y = rk4_step(f, x, y, dx)
    return x, y
    
def main():

    dx = 0.01
    x, y = integrate(f, 0.0, 0.0, 1.0, 2.0, dx)
    print 'Solution at x =', x, 'is', y

if __name__ in ('__main__'):
    main()
