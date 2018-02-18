import sys
import math
import numpy as np
import matplotlib.pyplot as plt

def initialize(n, seed, v0, rj, mj):
    mass = np.ones(n)/n
    pos = np.zeros((n,3))
    vel = np.zeros((n,3))
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
    elif n==3:#sun, earth, jupiter
        mass = np.array([1.0, 3e-6, mj])
        ve = math.sqrt((mass[0] + mass[1])/1.0)#re=1
        vj = math.sqrt((mass[0] + mass[2])/rj)

        pos = np.array([[0.,0,0],[1.0,0,0],[0,rj,0]])
        vel = np.array([[0,0,0],[0,ve,0],[-vj,0,0]])
    else:
        r = (np.random.random(n))**(1./3)
        theta = np.arccos(2*np.random.random(n)-1)
        phi = 2*np.pi*np.random.random(n)
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        pos = zip(x,y,z)

        v = 0.7
        theta = np.arccos(2*np.random.random(n)-1)
        phi = 2*np.pi*np.random.random(n)
        x = v*np.sin(theta)*np.cos(phi)
        y = v*np.sin(theta)*np.sin(phi)
        z = v*np.cos(theta)
        vel = zip(x,y,z)
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

def initial_pythag(seed, p2y):
    mass = np.array([5.0, 4.0, 3.0])
    pos = np.array([[0,0,0],[0,3.0,0],[4.0,p2y,0]])
    vel = np.array([[0,0,0],[0,0,0],[0,0,0.0]])
    return mass,pos,vel

def rms_dist(mass, pos1, pos2):
    temp = mass*((pos1-pos2)**2).sum(axis=1)
    M_inv = 1.0/(mass.sum())
    return math.sqrt(M_inv*temp.sum())

def rms_size(mass, pos):
    temp = mass*(pos**2).sum(axis=1)
    M_inv = 1.0/(mass.sum())
    return math.sqrt(M_inv*temp.sum())

def plot_snap(t, pos):
    plt.figure()
    plt.title("Time t = %s"%t)
    posT = np.transpose(pos)
    plt.plot(posT[0], posT[1], '.')
    plt.savefig("3c-time%s.png"%t)
    plt.clf()
