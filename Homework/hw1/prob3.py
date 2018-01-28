#Problem 3
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import math

global rp, rm
global E, GM, L, epsilon

def phi(r):
    return -GM/(r*r + epsilon**2)**0.5

def v_sqr(r):
    E = -0.5
    L = 0.5
    return (2*(E-phi(r))-(L/r)**2)

def f(r):
    return (2*(r-rm)**0.5 * (rp-r)**0.5)/(v_sqr(r))**0.5


def gauss_cheb(f, n):

    x = []
    for i in range(n):
        x.append(math.cos((2*i+1)*math.pi/(2.*n)))
    w = math.pi/n

    sum = 0.0
    for i in range(len(x)):
        r = (rp + rm)/2 + x[i]*(rp - rm)/2
        sum += w*f(r)
    return sum

#define consts
E = -0.5
L = 0.5
epsilon = 0.1
GM = 1

#find roots
rm = float(fsolve(v_sqr, -2))
rp = float(fsolve(v_sqr, 2))


P = gauss_cheb(f, 10)*(rp - rm)/2
print "epsilon = 0.1, period = ", P


#-------part b-------
#re-define epsilon
epsilon = 0

#find roots
rm = float(fsolve(v_sqr, -2))
rp = float(fsolve(v_sqr, 2))


P = gauss_cheb(f, 10)*(rp - rm)/2
print "epsilon = 0, period = ", P

newtonian = 2*math.pi*GM*(-2*E)**(-1.5)
print "Newtonian period = ", newtonian
print "diff = ", abs(P-newtonian)
