from functions import *

def outputS(t, E0, mass, pos, vel, eps2, steps):
    E = energy(mass, pos, vel, eps2)
    print 't =', t, 'dE =', E-E0, 'steps =', steps


def main(N, seed, eps, dt, t_end, v0, dEtol):
    if eps <= 0.0: eps = tiny
    #dt0 = dt

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, seed, v0)
    steps = 0

    # Initial diagnostics.

    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    outputS(t, E0, mass, pos, vel, eps**2, steps)
    a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0], pos[1],
                                  vel[0], vel[1], eps**2)
    print 'semimajor axis =', a, ' eccentricity =', e

    # Run forward to specified time.

    tplot = []
    dEplot = []
    hplot = []
    smaplot = []
    eccplot = []

    while t < t_end-0.5*dt:
        
        t2 = t
        pos2 = pos.copy()
        vel2 = vel.copy()

        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        E1 = energy(mass, pos, vel, eps**2)

        t2,pos2,vel2 = step(t2, mass, pos2, vel2, eps**2, dt)
        E2 = energy(mass, pos2, vel2, eps**2)

        DeltaE = abs(E2 - E1)
        dtnext = dt*(DeltaE/dEtol)**(-1./3)

        if dtnext > 1.5*dt:
            dt = 1.5*dt
        else:
            dt = dtnext

        steps += 1
        E = energy(mass, pos, vel, eps**2)
        a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0],
                                      pos[1], vel[0], vel[1], eps**2)

        tplot.append(t)
        dEplot.append(abs(E-E0))
        hplot.append(h)
        smaplot.append(a)
        eccplot.append(e)
        
        if (abs(E-E0) > 1e-4):
            print "Energy error too great, breaking loop"
            print "dE = ", abs(E-E0)
            break
    # Final diagnostics.

    outputS(t, E0, mass, pos, vel, eps**2, steps)
    print "dEtol = ", dEtol
    print "Max |E-E0| = ", max(dEplot)

    plt.figure()

    plt.subplot(2,2,1)
    plt.plot(tplot, dEplot)
    plt.xlabel('time')
    plt.ylabel('energy error')

    plt.subplot(2,2,2)
    plt.plot(tplot, hplot)
    plt.xlabel('time')
    plt.ylabel('angular momentum')

    plt.subplot(2,2,3)
    plt.plot(tplot, smaplot)
    plt.xlabel('time')
    plt.ylabel('semimajor axis')

    plt.subplot(2,2,4)
    plt.plot(tplot, eccplot)
    plt.xlabel('time')
    plt.ylabel('eccentricity')

    plt.tight_layout()
    plt.savefig("2.png")
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
                      dest="eps", type="float", default ="0.0",
                      help="softening length eps [%default]")
    result.add_option("-d",
                      dest="dt", type="float", default ="0.01",
                      help="time step [%default]")
    result.add_option("-t",
                      dest="t_end", type="float", default ="50.0",
                      help="integration interval [%default]")
    result.add_option("-v",
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    result.add_option("-q",
                    dest="dEtol", type="float", default ="2e-6",
                    help="difference in energy [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
