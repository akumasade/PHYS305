from functions import *

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0, dEtol):
    if eps <= 0.0: eps = tiny
    np.random.seed(12345)
    # Initial conditions.
    t = 0.0
    steps = 0

    #cluster 1
    mass1,pos1,vel1 = initialize(N, seed, v0)
    pos1 = pos1 + np.array([2, 0.5, 0])
    vel1 = vel1 + np.array([-0.5, 0.0,0])

    #cluster 2
    mass2,pos2,vel2 = initialize(N, seed, v0)
    pos2 = pos2 + np.array([-2, -0.5, 0])
    vel2 = vel2 + np.array([0.5, 0.0,0])

    #combinte them
    mass = np.concatenate((mass1, mass2), axis=0)
    pos = np.concatenate((pos1, pos2), axis=0)
    vel = np.concatenate((vel1, vel2), axis=0)


    # Initial diagnostics.
    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2)

    # Run forward to specified time.
    tplot = []
    dEplot = []

    Eplot = []
    Uplot = []
    Rplot = []

    time_interval = np.array([x for x in range(int(math.ceil(t_end+3)))])
    ind = 1

    while t < t_end-0.5*dt:
        #apply variable time step (step doubling)

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
        #a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0],
        #                              pos[1], vel[0], vel[1], eps**2)

        U = potential_energy(mass, vel, eps**2)
        R = virial_rad(mass, pos, eps**2)
        tplot.append(t)
        dEplot.append(abs(E-E0))

        Eplot.append(E)
        Uplot.append(U)
        Rplot.append(R)
        #print t
        #print energy each time unit
        if t>=time_interval[ind]-dt and t<=time_interval[ind]+dt:
            print "t=", t
            print "dE = ",abs(E-E0)
            ind +=1
            
        if (abs(E-E0) > 1e-4):
            print "Energy error too great, breaking loop"
            print "dE = ", abs(E-E0)
            break

    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2)
    #print "Max eccentricity:", max(eccplot)
    print "dEtol = ", dEtol
    print "Max |E-E0| = ", max(dEplot)

    plt.plot(tplot, Eplot, label="Total Energy")
    plt.plot(tplot, Uplot, label="Potential Energy")
    plt.plot(tplot, Rplot, label="Virial Radius")
    plt.xlabel('t')
    plt.ylabel('E')
    plt.legend()
    plt.savefig("3.png")
    plt.clf()


def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-n",
                      dest="N", type="int", default ="50",#10 for now, but goddamn this is gonna suck
                      help="number of particles [%default]")
    result.add_option("-s",
                      dest="seed", type="int", default ="12345",
                      help="random seed [%default]")
    result.add_option("-e",
                      dest="eps", type="float", default ="0.001",
                      help="softening length eps [%default]")
    result.add_option("-d",
                      dest="dt", type="float", default ="0.001",
                      help="time step [%default]")
    result.add_option("-t",
                      dest="t_end", type="float", default ="30.0",
                      help="integration interval [%default]")
    result.add_option("-v",
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    result.add_option("-q",
                    dest="dEtol", type="float", default ="5e-8",
                    help="difference in energy [%default]")

    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
