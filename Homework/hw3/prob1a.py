from functions import *

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0, rj, mj):
    if eps <= 0.0: eps = tiny

    # Initial conditions.

    t = 0.0
    mass,pos,vel = initialize(N, seed, v0, rj, mj)

    # Initial diagnostics.

    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2)
    a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0], pos[1],
                                  vel[0], vel[1], eps**2)
    print 'semimajor axis =', a, ' eccentricity =', e

    # Run forward to specified time.

    tplot = []
    dEplot = []
    hplot = []
    smaplot = []
    eccplot = []

    rp = 1.e6
    rpp = rp+1.
    posp = np.zeros(3)
    pospp = np.ones(3)
    while t < t_end-0.5*dt:
        t,pos,vel = step(t, mass, pos, vel, eps**2, dt)
        E = energy(mass, pos, vel, eps**2)
        a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0],
                                      pos[1], vel[0], vel[1], eps**2)
        r = (((pos[0]-pos[1])**2).sum())**0.5
        #print t, E-E0, a, e, r
        if r < rp and rp >= rpp:
            v1 = (rp-rpp)/dt
            v2 = (r-rp)/dt
            tmax = t - 1.5*dt + dt*(-v1)/(v2-v1)
            xmax = 0.5*(pospp+posp) + 0.5*(pos[0]-pos[1]-pospp)*(-v1)/(v2-v1)
            print 'maximum', tmax, rp, math.atan2(posp[1], posp[0]), \
                  (xmax**2).sum()**0.5, math.atan2(xmax[1], xmax[0]) \

        tplot.append(t)
        dEplot.append(E-E0)
        hplot.append(h)
        smaplot.append(a)
        eccplot.append(e)

        rpp = rp
        rp = r
        pospp = posp.copy()
        posp = pos[0]-pos[1]

    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2)
    print "Max eccentricity:", max(eccplot)
    #only need gto plot eccentricity
    plt.figure()

    #plt.subplot(2,2,4)
    plt.plot(tplot, eccplot)
    plt.xlabel('time')
    plt.ylabel('eccentricity')
    plt.title('$M_{j}=$%s, $R_{j}=$%s'%(mj, rj))

    #plt.tight_layout()
    #plt.show()
    plt.savefig("1a.png")

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
