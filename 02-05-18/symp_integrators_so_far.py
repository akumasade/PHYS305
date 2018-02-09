# Symplectic integrators are specific to Hamiltonian systems, where
# the acceleration acc is a function of position only.  These examples
# are written using "dynamical" notation, so the independent variable
# is time t, pos is position, and vel is velocity.

# Second order predictor-corrector.

def pc2_step(acc, t, pos, vel, dt):
    a0 = acc(pos)
    pos += (vel + 0.5*a0*dt)*dt
    a1 = acc(pos)
    vel += 0.5*(a0+a1)*dt
    t += dt
    return t, pos, vel

# Kick-drift-kick, 2nd-order.

def kdk_step(acc, t, pos, vel, dt):
    vel += 0.5*acc(pos)*dt
    pos += vel*dt
    vel += 0.5*acc(pos)*dt
    t += dt
    return t, pos, vel

# Drift-kick-drift, 2nd-order.

def dkd_step(acc, t, pos, vel, dt):
    pos += 0.5*vel*dt
    vel += acc(pos)*dt
    pos += 0.5*vel*dt
    t += dt
    return t, pos, vel
