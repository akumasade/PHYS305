# nbody_ex5.py: starting point for a simple N-body code.
# example 6.2
# Initializes, implements acceleration/potential calculation,
# prints a simple diagnostic message, and runs to completion.

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def initialize(n, seed, v0, rj, mj):
    mass = np.ones(n)/n
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
    elif n==3:#sun, earth, jupiter
        mass = np.array([1.0, 3e-6, mj])
        ve = math.sqrt((mass[0] + mass[1])/1.0)#re=1
        vj = math.sqrt((mass[0] + mass[2])/rj)

        pos = np.array([[0.,0,0],[1.0,0,0],[rj,0,0]])
        vel = np.array([[0,0,0],[0,ve,0],[-vj,0,0]])
    else:
        np.random.seed(12345)
        r = (np.random.random(n))**(1./3)
        theta = np.arccos(2*np.random.random(n)-1)
        phi = 2*np.pi*np.random.random(n)
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        pos = zip(x,y,z)

        vel = np.zeros((n,3))
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

def kinetic_energy(mass, vel):
    return 0.5*(mass*(vel**2).sum(axis=1)).sum()

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


def step(t, mass, pos, vel, eps2, dt):

    # Second-order predictor-corrector.

    acc = acceleration(mass, pos, eps2)
    pos += dt*(vel+0.5*dt*acc)
    anew = acceleration(mass, pos, eps2)
    vel += 0.5*dt*(acc+anew)

    return t+dt,pos,vel

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
def main(N, seed, eps, dt, t_end, v0, mj, rj):
    if eps <= 0.0: eps = tiny

    # Initial conditions.
    t = 0.0
    mass,pos,vel = initialize(N, seed, v0, rj, mj)

    # Initial diagnostics.

    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2)

    #Plotting
    Es=[]
    ts = []
    saxs = []
    es = []
    rprev = pos[0]-pos[1]
    M = mass.sum()


    # Run forward to specified time.
    while t < (t_end-0.5*dt):
        r = pos[0]-pos[1]
        v = vel[0]-vel[1]
        h = np.linalg.norm(np.cross(r,v))
        r_norm = np.linalg.norm(r)

        E = energy(mass, pos, vel, eps**2)
        earth = orbital_elements(mass[0], mass[1], pos[0], pos[1], vel[0], vel[1], eps**2)
        #returns sma, ecc, E, h2**2

        #for plotting
        Es.append(E-E0)
        ts.append(t)
        saxs.append(earth[0])
        es.append(earth[1])
        if earth[0]<=0: break

        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)

        #finding maxima
        '''
        rnext = np.linalg.norm(pos[0]-pos[1])
        if (rnext-r_norm<=0) and (r_norm-np.linalg.norm(rprev)>=0):
            print "t = ",t, "max = ", r
            print "angle = ", np.arctan(r[1]/r[0])
        rprev = r
        '''

    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2)
    print "Max ecc:", max(es)


    fig = plt.figure()
    subplot1 = fig.add_subplot(2, 1, 1)
    plt.plot(ts, es)
    #plt.legend()
    plt.xlabel('t')
    plt.ylabel("eccentricity")

    subplot1 = fig.add_subplot(2, 1, 2)
    #plt.plot(ts, es, label="eccentricity")
    plt.plot(ts, saxs)
    #plt.legend()
    plt.xlabel('t')
    plt.ylabel("semi-major axis")
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
                      dest="eps", type="float", default ="0.00",
                      help="softening length eps [%default]")
    result.add_option("-d",
                      dest="dt", type="float", default ="0.01",
                      help="time step [%default]")
    result.add_option("-t",
                      dest="t_end", type="float", default ="1000.0",
                      help="integration interval [%default]")
    result.add_option("-v",
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    result.add_option("-r",
                      dest="rj", type="float", default ="2.5",
                      help="initial jubpiter radius [%default]")
    result.add_option("-m",
                      dest="mj", type="float", default ="0.05",
                      help="jupiter mass [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
