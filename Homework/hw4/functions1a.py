import sys
import math
import numpy as np
import matplotlib.pyplot as plt

tiny = 1.e-20
huge = 1.e20

def initialize(n, seed, v0):
    np.random.seed(seed)
    mass = np.ones(n)/n
    pos = np.zeros((n,3))
    vel = np.zeros((n,3))
    if n == 2:
        pos = np.array([[1.0,0,0],[-1.0,0,0]])
        vel = np.array([[0,v0,0],[0,-v0,0]])
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

'''
def potential_energy(mass, pos, eps2):
    n = len(mass)
    pot = 0.0
    dx = np.zeros((n,3))
    for i in range(n):
        dx[i+1:] = pos[i+1:] - pos[i]
        dr2 = (dx[i+1:]**2).sum(axis=1) + eps2
        pot -= mass[i]*(mass[i+1:]/np.sqrt(dr2)).sum()
    return pot
'''
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

def output(t, E0, mass, pos, vel, eps2, steps):
    E = energy(mass, pos, vel, eps2)
    print 't =', t, 'dE =', E-E0, 'steps =', steps

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

def step(t, mass, pos, vel, eps2, dt):

    # Second-order predictor-corrector.

    acc,tau = acceleration2(mass, pos, vel, eps2)
    pos += dt*(vel+0.5*dt*acc)
    anew,tau = acceleration2(mass, pos, vel, eps2)
    vel += 0.5*dt*(acc+anew)
    return t+dt,pos,vel,tau

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

def virial_rad(mass, pos, eps2):
    U = potential_energy(mass, pos, eps2)
    return -(mass.sum())**2 / (2*U)
