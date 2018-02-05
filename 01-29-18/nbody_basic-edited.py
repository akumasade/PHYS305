# nbody_basic.py: starting point for a simple N-body code.
# Runs, but does nothing...

import sys
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, seed):
    np.random.seed(seed)
    v0 = 0.25
    x0 = 1.0
    mass = np.ones(n)/n#masses sum up to 1
    if n==2:
        pos = np.array([[1,0,0], [-1,0,0]])*x0
        vel = np.array([[0,1,0], [0,-1,0]])*v0
    else:
        pos = np.zeros((n,3))
        vel = np.zeros((n,3))

    return mass,pos,vel

def get_pot(mass, pos, eps2):
    n = len(mass)

    U = 0
    dx = np.zeros((n,3))
    for i in range(n):
        dx[i+1:,:] = pos[i+1:,:] - pos[i,:]
        dr2 = (dx[i+1:,:]**2).sum(axis=1) + eps2
        U -= mass[i]*(mass[i+1:]/np.sqrt(dr2)).sum()
    #print U
    '''
    pot = 0
    for i, m in enumerate(mass):
        pot += -m*mass[i:]/np.sqrt((pos[i,:]-pos[i:,:])**2 + eps2*eps2)
    print pot
    '''
    return U

def energy(mass, pos, vel, eps2):
    T = np.sum(0.5*mass*(vel**2).sum())#prob wrong but whatever
    U = get_pot(mass, pos, eps2)
    return T+U

def output(t, E0, mass, pos, vel, eps2):
    E = energy(mass, pos, vel, eps2)
    print 't =', t, 'dE =', E-E0

def get_acc(mass, pos, eps2):
    n = len(mass)
    acc = np.zeros((n,3))
    dx = np.zeros((n,3))
    for i in range(n):
        dx[i+1:,:] = pos[i+1:,:] - pos[i,:]
        dr2 = (dx[i+1:,:]**2).sum(axis=1)**2
        acc -= mass[i]*(mass[i+1:]/np.sqrt(dr2 +eps2*eps2))

    return acc

def step(t, mass, pos, vel, eps2, dt):
    return t+dt,pos,vel

def main(N, seed, eps2, dt, t_end):
    eps2 = eps2**2

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, seed)

    # Initial diagnostics.
    E0 = energy(mass, pos, vel, eps2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps2)

    # Run forward to specified time.

    while t < t_end-0.5*dt:
        t,pos,vel = step(t, mass, pos, vel, eps2, dt)

    # Final diagnostics.
    output(t, E0, mass, pos, vel, eps2)

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
                      dest="eps2", type="float", default ="0.01",
                      help="softening length eps [%default]")
    result.add_option("-d",
                      dest="dt", type="float", default ="0.01",
                      help="time step [%default]")
    result.add_option("-t",
                      dest="t_end", type="float", default ="10.0",
                      help="integration interval [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
