
def g(z):				# algebraic equation to solve
    return z*z - 2.

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

def main():
    zl = 1.0
    zr = 2.0
    n,zl,zr = bisect(g, zl, zr, 1.e-6)

    print "Root lies in range (%f, %f) after %d iterations"%(zl, zr, n)
    print "Function value =", g(0.5*(zl+zr))

if __name__ in ('__main__'):
    main()
