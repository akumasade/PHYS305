# Solution to Exercise 6.1.

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, seed, v0):
    mass = np.ones(n)/n
    pos = np.zeros((n,3))
    vel = np.zeros((n,3))
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
    return mass,pos,vel

def potential_energy(mass, pos, eps2):
    n = len(mass)
    pot = 0.0
    dx = np.zeros((n,3))
    for i in range(n):
        dx[i+1:] = pos[i+1:] - pos[i]
        dr2 = (dx[i+1:]**2).sum(axis=1) + eps2
        pot -= mass[i]*(mass[i+1:]/np.sqrt(dr2)).sum()
    return pot

def alt_potential_energy(mass, pos, eps2):
    n = len(mass)
    pot = 0.0
    for i in range(n):
        for j in range(i+1,n):
            dr2 = eps2
            for k in range(3):
                dr2 += (pos[i,k]-pos[j,k])**2
            pot -= mass[i]*mass[j]/math.sqrt(dr2)
    return pot

def kinetic_energy(mass, vel):
    return 0.5*(mass*(vel**2).sum(axis=1)).sum()

def alt_kinetic_energy(mass, vel):
    n = len(mass)
    kin = 0.0
    for i in range(n):
        vi2 = 0.0
        for k in range(3):
            vi2 += vel[i,k]**2
        kin += 0.5*mass[i]*vi2
    return kin

def energy(mass, pos, vel, eps2):
    T = kinetic_energy(mass, vel)
    U = potential_energy(mass, pos, eps2)
    return T+U

def output(t, E0, mass, pos, vel, eps2):
    E = energy(mass, pos, vel, eps2)
    print 't =', t, 'dE =', E-E0

def acceleration(mass, pos, eps2):
    n = len(mass)
    acc = np.zeros((n,3))
    for i in range(n):
        dx   = pos - pos[i]
        dr2  = (dx**2).sum(axis=1) + eps2
        dr2i = 1./dr2
        dr3i = mass*np.sqrt(dr2i)*dr2i
        dx  *= dr3i.reshape(n,1)
        acc[i] = dx.sum(axis=0)
    return acc

def alt_acceleration(mass, pos, eps2):
    n = len(mass)
    acc = np.zeros((n,3))
    for i in range(n):
        for j in range(i+1,n):
            dr2 = eps2
            for k in range(3):
                dr2 += (pos[j,k]-pos[i,k])**2
            dr2i = 1./dr2
            dr3i = dr2i*math.sqrt(dr2i)
            for k in range(3):
                dxij = (pos[j,k]-pos[i,k])*dr3i
                acc[i,k] += mass[j]*dxij
                acc[j,k] -= mass[i]*dxij
    return acc

def step(t, mass, pos, vel, eps2, dt):

    # Second-order predictor-corrector.

    acc = acceleration(mass, pos, eps2)
    pos += dt*(vel+0.5*dt*acc)
    anew = acceleration(mass, pos, eps2)
    vel += 0.5*dt*(acc+anew)

    return t+dt,pos,vel

def kdk_step(t, mass, pos, vel, eps2, dt):
    vel += 0.5*acceleration(mass, pos, eps2)*dt
    pos += vel*dt
    vel += 0.5*acceleration(mass, pos, eps2)*dt
    t += dt
    return t, pos, vel

def orbital_elements(m1, m2, x1, x2, v1, v2, eps2):
    M = m1+m2
    x = x2-x1
    v = v2-v1
    r2 = (x**2).sum() + eps2
    v2 = (v**2).sum()
    E = 0.5*v2 - M/math.sqrt(r2)
    sma = -0.5*M/E
    h2 = ((np.cross(x,v))**2).sum()
    ecc = (1 + 2*E*h2/M**2)**0.5
    return sma, ecc, E, h2**0.5
    
tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0):
    if eps <= 0.0: eps = tiny

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, seed, v0)

    # Initial diagnostics.
    
    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2)
    a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0], pos[1],
                                  vel[0], vel[1], eps**2)
    print 'semimajor axis =', a, ' eccentricity =', e

    # Run forward to specified time.

    tplot = []
    dEplot = []
    hplot = []
    smaplot = []
    eccplot = []

    rp = 1.e6
    rpp = rp+1.
    posp = np.zeros(3)
    pospp = np.ones(3)
    while t < t_end-0.5*dt:
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        E = energy(mass, pos, vel, eps**2)
        a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0],
                                      pos[1], vel[0], vel[1], eps**2)
        r = (((pos[0]-pos[1])**2).sum())**0.5
        #print t, E-E0, a, e, r
        if r < rp and rp >= rpp:
            v1 = (rp-rpp)/dt
            v2 = (r-rp)/dt
            tmax = t - 1.5*dt + dt*(-v1)/(v2-v1)
            xmax = 0.5*(pospp+posp) + 0.5*(pos[0]-pos[1]-pospp)*(-v1)/(v2-v1)
            print 'maximum', tmax, rp, math.atan2(posp[1], posp[0]), \
                  (xmax**2).sum()**0.5, math.atan2(xmax[1], xmax[0]) \

        tplot.append(t)
        dEplot.append(E-E0)
        hplot.append(h)
        smaplot.append(a)
        eccplot.append(e)

        rpp = rp
        rp = r
        pospp = posp.copy()
        posp = pos[0]-pos[1]

    # Final diagnostics.
    
    output(t, E0, mass, pos, vel, eps**2)

    plt.figure()

    plt.subplot(2,2,1)
    plt.plot(tplot, dEplot)
    plt.xlabel('time')
    plt.ylabel('energy error')

    plt.subplot(2,2,2)
    plt.plot(tplot, hplot)
    plt.xlabel('time')
    plt.ylabel('angular momentum')

    plt.subplot(2,2,3)
    plt.plot(tplot, smaplot)
    plt.xlabel('time')
    plt.ylabel('semimajor axis')

    plt.subplot(2,2,4)
    plt.plot(tplot, eccplot)
    plt.xlabel('time')
    plt.ylabel('eccentricity')

    plt.tight_layout()
    plt.show()

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", 
                      dest="N", type="int", default ="2",
                      help="number of particles [%default]")
    result.add_option("-s", 
                      dest="seed", type="int", default ="42",
                      help="random seed [%default]")
    result.add_option("-e", 
                      dest="eps", type="float", default ="0.0",
                      help="softening length eps [%default]")
    result.add_option("-d", 
                      dest="dt", type="float", default ="0.01",
                      help="time step [%default]")
    result.add_option("-t", 
                      dest="t_end", type="float", default ="50.0",
                      help="integration interval [%default]")
    result.add_option("-v", 
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    return result
    
if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
