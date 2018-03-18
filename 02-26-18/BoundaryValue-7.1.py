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

def g(z):		# algebraic equation to solve
    x, y = integrate(f, 0, 1.0, z, 2.0, 0.01 )
    return y - y1

def bisect(func, zl, zr, tol):
    n = 0
    while zr-zl > tol:
        zm = 0.5*(zl + zr)
    	if func(zm)*func(zl) > 0:
    	    zl = zm
    	else:
    	    zr = zm
    	n += 1
        print n, zl, func(zl), zr, func(zr)

    return n,zl,zr

def falsepos(z1,z2):
    zs = z1 + (z2 - z1) * (-g(z1)) / (g(z2) - g(z1))

    if (g(zs)*g(z1) > 0):
        z1 = zs
    else:
        z2 = zs
    return zs

def main():
    global y0, y1
    y0 = 1.0
    y1 = 2.0
    y = [y0, y1]#y vals
    x = [0.0, 1.0]#corresponding times
    #initial z guess
    z0 = 3
    g0 = g(z0)

    z1 = -1
    g1 = g(z1)

    print g0, g1
    #perform bisection
    tol = 1e-4#tolerance
    low = min([g1, g0])
    high = max([g1, g0])
    z = bisect(f, low, high, tol)

    print z
    #perform false-position


    #correct answer




if __name__ in ('__main__'):
    main()
