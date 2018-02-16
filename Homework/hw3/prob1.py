from functions import *
import matplotlib.animation as animation

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0, rj, mj):
    if eps <= 0.0: eps = tiny
    global t, mass, pos, vel, t_dia, t_anim
    # Initial conditions.
    t = 0.0
    mass,pos,vel = initialize(N, seed, v0, rj, mj)
    print vel

    # Move to center of mass frame
    m = mass.reshape(N,1)/mass.sum()
    pos -= (m*pos).sum(axis=0)
    vel -= (m*vel).sum(axis=0)

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

    '''
    # Run forward to specified time.
    while t < (t_end-0.5*dt):

        E = energy(mass, pos, vel, eps**2)
        earth = orbital_elements(mass[1], mass[0], pos[1], pos[0], vel[1], vel[0], eps**2)
        #returns sma, ecc, E, h2**2

        #for plotting
        Es.append(E-E0)
        ts.append(t)
        saxs.append(earth[0])
        es.append(earth[1])
        if earth[0]<=0: break

        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)

        if t-np.floor(t)==0: print t
    '''

    dt_anim = 0.05
    dt_dia = 1.0
    t = 0.0
    t_dia = 0.0
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
            E = energy(mass, pos, vel, eps**2)
            earth = orbital_elements(mass[1], mass[0], pos[1], pos[0], vel[1], vel[0], eps**2)
            #returns sma, ecc, E, h2**2

            #for plotting
            Es.append(E-E0)
            ts.append(t)
            saxs.append(earth[0])
            es.append(earth[1])

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

    # Final diagnostics.

    #output(t, E0, mass, pos, vel, eps**2)
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
                      dest="N", type="int", default ="3",
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
