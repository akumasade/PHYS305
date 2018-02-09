# nbody_ex3.py: starting point for a simple N-body code.

# Initializes, implements acceleration/potential calculation, and
# prints a simple diagnostic message.

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, seed, v0):
    mass = np.ones(n)/n
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
    else:
        pos = np.zeros((n,3))
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
    return t+dt,pos,vel

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

    acc = acceleration(mass, pos, eps**2)
    print 'initial acc =', acc
    acc = alt_acceleration(mass, pos, eps**2)
    print 'initial alt_acc =', acc

    '''
    # Run forward to specified time.

    while t < t_end-0.5*dt:
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)

    # Final diagnostics.
    
    output(t, E0, mass, pos, vel, eps**2)

    plt.figure(figsize=(6,6))
    plt.scatter(pos[:,0], pos[:,1])
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Exercise 3')
    plt.show()
    '''

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
                      dest="eps", type="float", default ="0.05",
                      help="softening length eps [%default]")
    result.add_option("-d", 
                      dest="dt", type="float", default ="0.001",
                      help="time step [%default]")
    result.add_option("-t", 
                      dest="t_end", type="float", default ="10.0",
                      help="integration interval [%default]")
    result.add_option("-v", 
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    return result
    
if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
