from functions import *

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0, rj, mj):
    if eps <= 0.0: eps = tiny
    np.random.seed(12345)
    # Initial conditions.
    t = 0.0
    time_interval = np.array([1.0*x for x in range(int(math.ceil(t_end+1)))])
    ind = 1

    #cluster 1
    mass1,pos1,vel1 = initialize(N, seed, v0, rj, mj)
    pos1 = pos1 + np.array([2, 0.5, 0])
    vel1 = vel1 + np.array([-0.5, 0.0,0])
    mass1 = mass1*0.5
    #cluster 2
    mass2,pos2,vel2 = initialize(N, seed, v0, rj, mj)
    pos2 = pos2 + np.array([-2, -0.5, 0])
    vel2 = vel2 + np.array([0.5, 0.0,0])
    mass2 = mass2*0.5

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
    KEplot = []
    RMSplot = []


    while t < t_end-0.5*dt:

        #step doubling
        t0 = t
        pos0 = pos.copy()
        vel0 = vel.copy()

        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        E1 = energy(mass, pos, vel, eps**2)

        t2,pos2,vel2 = step(t0, mass, pos0, vel0, eps**2, dt)
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

        KE = kinetic_energy(mass, vel)
        RMS = rms_size(mass, pos)
        tplot.append(t)
        dEplot.append(E-E0)

        Eplot.append(E)
        KEplot.append(KE)
        RMSplot.append(RMS)

        if t>time_interval[ind]-dt and t<=time_interval[ind]:
            print "t = ", t
            print "|E-E0|",abs(E-E0)
            ind +=1


    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2)
    #print "Max eccentricity:", max(eccplot)


    plt.plot(tplot, Eplot)
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title('Total Energy plot')
    plt.legend()
    plt.savefig("3c-E.png")
    plt.clf()

    plt.plot(tplot, KEplot, 'r')
    plt.xlabel('t')
    plt.ylabel('KE')
    plt.title('Kinetic Energy plot')
    plt.legend()

    plt.savefig("3c-KE.png")
    plt.clf()

    plt.plot(tplot, RMSplot)
    plt.xlabel('t')
    plt.ylabel('RMS')
    plt.title('RMS Size')

    #plt.tight_layout()
    #plt.show()
    plt.savefig("3c-RMS.png")
    plt.clf()

def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-n",
                      dest="N", type="int", default ="150",#10 for now, but goddamn this is gonna suck
                      help="number of particles [%default]")
    result.add_option("-s",
                      dest="seed", type="int", default ="12345",
                      help="random seed [%default]")
    result.add_option("-e",
                      dest="eps", type="float", default ="0.1",
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
    result.add_option("-r",
                      dest="rj", type="float", default ="2.5",
                      help="initial jubpiter radius [%default]")
    result.add_option("-m",
                      dest="mj", type="float", default ="0.05",
                      help="jupiter mass [%default]")
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
