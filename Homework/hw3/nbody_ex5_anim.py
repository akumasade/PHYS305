# nbody_ex5.py: starting point for a simple N-body code.

# Initializes, implements acceleration/potential calculation, 
# prints a simple diagnostic message, and runs to completion.

import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def initialize(n, seed, v0, eps2):
    mass = np.ones(n)/n
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
    elif n == 3:

        # Add a third body.
            
        pos = np.array([[1.0,0,0],[-1.0,0,0],[0,0.5,0]])
        vel = np.array([[0,v0,0],[0,-v0,0],[0.,0,0]])

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

tiny = 1.e-20
def main(N, seed, eps, dt, dt_dia, dt_anim, t_end, v0):
    global t, mass, pos, vel, t_dia, t_anim
    if eps <= 0.0: eps = tiny

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, seed, v0, eps**2)

    # Initial diagnostics.
    
    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2)
    t_dia = dt_dia

    # Initialize animation.
    
    fig = plt.figure(figsize=(6,6))
    scat = plt.scatter(pos[:,0], pos[:,1], s=20, c=['r','g','b'])
    lim = 1.25*max(abs(pos[:,0].min()), abs(pos[:,0].max()),
                   abs(pos[:,1].min()), abs(pos[:,1].max()))
    lim = 0.5*(int(lim/0.5)+1)
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.xlabel('x')
    plt.ylabel('y')
    text = plt.title('')
    text.set_text('time = {:7.1f}'.format(t))
    t_anim = dt_anim

    def update(frame):
        global t, mass, pos, vel, t_dia, t_anim

        # Run forward to specified time, under control of FuncAnimation.

        if t < t_end-0.5*dt:
            t,pos,vel = step(t, mass, pos, vel, eps**2, dt)

        if t >= t_dia-0.5*dt:
            t_dia += dt_dia
            output(t, E0, mass, pos, vel, eps**2)

        if t >= t_anim-0.5*dt:
            t_anim += dt_anim
            off = []
            for j in range(pos.shape[0]):
                off.append([pos[j,0], pos[j,1]])
                scat.set_offsets(off)
                text.set_text('time = {:7.1f}'.format(t))

        if t >= t_end-0.5*dt:
            anim.event_source.stop()

        return scat,

    anim = animation.FuncAnimation(fig, update, interval=1)
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
                      dest="eps", type="float", default ="0.05",
                      help="softening length eps [%default]")
    result.add_option("-d", 
                      dest="dt", type="float", default ="0.01",
                      help="time step [%default]")
    result.add_option("-A", 
                      dest="dt_anim", type="float", default ="0.05",
                      help="animation time step [%default]")
    result.add_option("-D", 
                      dest="dt_dia", type="float", default ="1.0",
                      help="diagnostic time step [%default]")
    result.add_option("-t", 
                      dest="t_end", type="float", default ="20.0",
                      help="integration interval [%default]")
    result.add_option("-v", 
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    return result
    
if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)
