# Solution to Homework 4, problem 1.
# Calculate time to reach stable KE
# (2 clusters of 50 particles)

import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, dx, v0, seed):
    mass = np.ones(n)/n
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])

    elif n == 3:

        # Pythagorean problem:
        mass = np.array([5.0, 4.0, 3.0])
        pos = np.array([[0.0,0,0],[1.0,0.000+dx,0],[0,1.0,0]])
        vel = np.array([[0.0,0,0],[0.0,0,0],[0.0,0,0]])

    else:

        np.random.seed(seed)
        rr = np.random.random(n)**(1./3)
        theta = np.arccos(2*np.random.random(n) - 1.)
        phi = 2*math.pi*np.random.random(n)
        x = rr*np.sin(theta)*np.cos(phi)
        y = rr*np.sin(theta)*np.sin(phi)
        z = rr*np.cos(theta)
        pos = np.array(zip(x, y, z))

        vv = np.ones(n)*v0
        theta = np.arccos(2*np.random.random(n) - 1.)
        phi = 2*math.pi*np.random.random(n)
        vx = vv*np.sin(theta)*np.cos(phi)
        vy = vv*np.sin(theta)*np.sin(phi)
        vz = vv*np.cos(theta)
        vel = np.array(zip(vx, vy, vz))

    # Move to the center of mass frame.
    
    m = mass.reshape(n,1)/mass.sum()
    pos -= (m*pos).sum(axis=0)
    vel -= (m*vel).sum(axis=0)

    if n >= 10:

        # Make a collision (doubles N, total mass).

        pos1 = pos.copy()
        vel1 = vel.copy()
        pos += np.array([1.,0.0,0])
        pos1 += np.array([-1.,0.0,0])
        vel += np.array([-0.5,0,0])
        vel1 += np.array([0.5,0,0])

        mass = np.concatenate((mass, mass))
        pos = np.concatenate((pos, pos1))
        vel = np.concatenate((vel, vel1))

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

huge = 1.e20
def acceleration2(mass, pos, vel, eps2):
    n = len(mass)
    acc = np.zeros((n,3))
    tau = huge
    for i in range(n):
        dx = pos - pos[i]
        dr2 = (dx**2).sum(axis=1) + eps2
        dv = vel - vel[i]
        vdotx = np.abs((dx*dv).sum(axis=1)) + tiny
        dr3 = dr2*np.sqrt(dr2)
        dr3i = mass/dr3
        dx *= dr3i.reshape(n,1)
        acc[i] = dx.sum(axis=0)

        dr3[i] = huge
        tau2 = dr2/vdotx
        tau3 = dr3/(mass+mass[i])
        tau2[i] = huge
        tau3[i] = huge
        tau = min(tau, tau2.min(), math.sqrt(tau3.min()))
    return acc, tau

def step2(t, mass, pos, vel, eps2, dt):
    '''
    # Second-order predictor-corrector.
    acc,tau = acceleration2(mass, pos, vel, eps2)
    pos += dt*(vel+0.5*dt*acc)
    anew,tau = acceleration2(mass, pos, vel, eps2)
    vel += 0.5*dt*(acc+anew)
    return t+dt,pos,vel,tau
    '''
    w = np.array([-1.17767998417887, 0.235573213359357, 0.784513610477560])
    wsum = 1.0 - 2.0*np.sum(w)
    
    c = [0.5*w[2], 0.5*(w[2] + w[1]), 0.5*(w[1] + w[0]), 0.5*(w[0] + wsum)]
    c = np.array(c + list(reversed(c))) #second half of c looks like fist half in reverse
    
    d = np.array([0.0, w[2], w[1], w[0], wsum, w[0], w[1], w[2]])
    
    for cc,dd in zip(dt*c, dt*d):
        vel += dd*acceleration2(mass, pos, vel, eps2)[0]
        pos += cc*vel
    
    tau = acceleration2(mass,pos, vel, eps2)[1]
    t += dt
    return t, pos, vel, tau

def print_bound_pairs(Escale, mass, pos, vel, eps2):
    n = len(mass)
    dx = np.zeros((n,3))
    dv = np.zeros((n,3))
    for i in range(n):
        dx[i+1:] = pos[i+1:] - pos[i]
        dv[i+1:] = vel[i+1:] - vel[i]
        dr2 = (dx[i+1:]**2).sum(axis=1) + eps2
        dv2 = (dv[i+1:]**2).sum(axis=1)
        mtot = mass[i] + mass[i+1:]
        mu = mass[i]*mass[i+1:]/mtot
        erel = mu*(0.5*dv2 - mtot/np.sqrt(dr2))
        bound = np.where(erel < -Escale)[0]
        if len(bound) > 0:
            for j in bound:
                s = 'bound pair: i = {:d}, j = {:d}, ' \
                      + 'Erel = {:9.5f} = {:8.3f} Escale'
                print s.format(i, j, erel[j], erel[j]/Escale)
    return

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

def plot_snap(t, pos, txt):
    plt.clf()
    plt.scatter(pos[:,0], pos[:,1], s=15)
    lim = 3.0
    plt.xlim(-lim, lim)
    plt.ylim(-lim, lim)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('t ={:5.1f}'.format(t))
    print 'Saved plot in', 'HW4.1 t ={:5.1f}.png'.format(t)
    plt.savefig('HW4.1'+txt+' t ={:5.1f}.png'.format(t))

tiny = 1.e-20
def main(N, seed, eps, dt_param, fixed_dt, t_end, dt_dia, v0):
    if eps <= 0.0: eps = tiny

    if fixed_dt:
        dt = dt_param
        txt = 'b'
    else:
        txt = 'a'

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, 0.0, v0, seed)
    E0 = energy(mass, pos, vel, eps**2)
    Escale = 2*abs(E0)/(3.*N)
    steps = 0

    if not fixed_dt:
        acc,tau = acceleration2(mass, pos, vel, eps**2)
        dt = dt_param*tau
    
    t_dia = dt_dia

    tplot = []
    Eplot = []
    Kplot = []
    Rplot = []
    Rvplot = []

    tsnap = [0.0, 1.0, 2.0, 5.0, 10.0, 30.0, 50.0, 1000000.0]

    fig = plt.figure(figsize=(6,6))
    #plot_snap(t, pos, txt)
    isnap = 1

    dE0 = 0.

    # Run forward to specified time.

    while t < t_end-0.5*dt:
        t,pos,vel,tau = step2(t, mass, pos, vel, eps**2, dt)
        steps += 1

        if not fixed_dt:
            dtnext = dt_param*tau
            if dtnext <= 1.5*dt:
                dt = dtnext
            else:
                print 'limiting dt:', t, dtnext, dt, dtnext/dt
                dt = 1.5*dt

        if t > t_dia-0.5*dt:
            t_dia += dt_dia
            E = energy(mass, pos, vel, eps**2)
            dE = E-E0
            print 't = {:6.3f}, dE = {:10.3e}, ddE = {:10.3e}, steps = {:6d}' \
                  .format(t, dE, dE-dE0, steps)
            dE0 = dE
            #print_bound_pairs(Escale, mass, pos, vel, eps**2)
            sys.stdout.flush()

        tplot.append(t)
        K = kinetic_energy(mass, vel)
        U = potential_energy(mass, pos, eps**2)
        Eplot.append(K+U)
        Kplot.append(K)
        R = math.sqrt((mass*((pos**2).sum(axis=1))).sum()/mass.sum())
        Rplot.append(R)
        Rvplot.append(-mass.sum()**2/(2.*U))
        
        
        #slope of kinetic energy over time
        if len(Kplot) > 1:
            Kslope = (Kplot[-2] - Kplot[-1])/dt
            #print Kslope
        else:
            Kslope = 100
        if (abs(Kslope) < 1e-4):#condition for "stable" KE
            print "Kinetic energy slope = ", Kslope, ", t = ", t
            

    plt.close()
    plt.figure()
    plt.plot(tplot, Eplot, label='E')
    plt.plot(tplot, Kplot, label='K')
    plt.plot(tplot, Rplot, label='R')
    plt.plot(tplot, Rvplot, label='Rv')
    plt.xlabel('time')
    plt.ylabel('E, K, R')
    plt.legend(loc='best')
    plt.savefig('1c.png')
    #plt.show()

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", 
                      dest="N", type="int", default ="50",
                      help="number of particles [%default]")
    result.add_option("-s", 
                      dest="seed", type="int", default ="12345",
                      help="random seed [%default]")
    result.add_option("-e", 
                      dest="eps", type="float", default ="0.001",
                      help="softening length eps [%default]")
    result.add_option("-d", 
                      dest="dt_param", type="float", default ="0.10",
                      help="time step parameter [%default]")
    result.add_option("-D", 
                      dest="dt_dia", type="float", default ="0.5",
                      help="diagnostic time step [%default]")
    result.add_option("-t", 
                      dest="t_end", type="float", default ="20.0",
                      help="integration interval [%default]")
    result.add_option("-v", 
                      dest="v0", type="float", default ="0.7",
                      help="initial 2-body v [%default]")
    result.add_option("-F", 
                      dest="fixed_dt", default = False, action='store_true',
                      help="fixed time step  dt_param[%default]")
    return result
    
if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
