# nbody_basic.py: starting point for a simple N-body code.
# Runs, but does nothing...

import sys
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, seed):
    np.random.seed(seed)
    mass = np.zeros(n)
    pos = np.zeros((n,3))
    vel = np.zeros((n,3))
    return mass,pos,vel
    
def get_pot(mass, pos, eps2):
    n = len(mass)
    pot = np.zeros(n)
    return pot

def energy(mass, pos, vel, eps2):
    T = 0.0
    U = 0.0
    return T+U

def output(t, E0, mass, pos, vel, eps2):
    E = energy(mass, pos, vel, eps2)
    print 't =', t, 'dE =', E-E0

def get_acc(mass, pos, eps2):
    n = len(mass)
    acc = np.zeros((n,3))
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
    from amuse.units.optparse import OptionParser
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
