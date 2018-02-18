from functions import *

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0, rj, mj):
    if eps <= 0.0: eps = tiny
    np.random.seed(12345)
    # Initial conditions.
    t = 0.0

    #cluster 1
    mass1,pos1,vel1 = initialize(N, seed, v0, rj, mj)
    pos1 = pos1 + np.array([2, 0.5, 0])
    vel1 = vel1 + np.array([-0.5, 0.0,0])
    mass1 = mass1*0.5
    #cluster 2

    mass2,pos2,vel2 = initialize(N, seed, v0, rj, mj)
    pos2 = pos2 + np.array([-2, -0.5, 0])
    vel2 = vel2 + np.array([-0.5, 0.0,0])
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
    snapshots = [0.0, 1.0, 2.0, 5.0, 10.0, 30.0]

    #plot t=0
    print t
    plot_snap(t, pos)

    while t < t_end-0.5*dt:
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)

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

        #plotting given snapshots
        if t in snapshots:
            print t
            plot_snap(t, pos)
    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2)
    #print "Max eccentricity:", max(eccplot)


    plt.plot(tplot, Eplot, label="Total")
    plt.plot(tplot, KEplot, label="Kinetic")
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title('Energy plot')
    plt.legend()

    #plt.tight_layout()
    #plt.show()
    plt.savefig("3c-energy.png")
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
                      dest="N", type="int", default ="50",#10 for now, but goddamn this is gonna suck
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
