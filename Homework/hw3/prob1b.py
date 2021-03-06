from functions import *
import matplotlib.patches as mpatches

tiny = 1.e-20
def main(N, seed, eps, dt, t_end, v0):
    if eps <= 0.0: eps = tiny

    #Create thing to keep track of our crap
    Mj = [0.01*x for x in range(1,21)]
    Rj = [0.1*x for x in range(12,26)]

    ecc_max = np.zeros((len(Mj),len(Rj)))

    #loop thru mj and rj
    for i,mj in enumerate(Mj):
        for j,rj in enumerate(Rj):

            # Initial conditions.

            t = 0.0
            mass,pos,vel = initialize(N, seed, v0, rj, mj)
            print mass
            print vel
            # Initial diagnostics.

            E0 = energy(mass, pos, vel, eps**2)
            #print 'Initial E =', E0
            #output(t, E0, mass, pos, vel, eps**2)
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
                    #print 'maximum', tmax, rp, math.atan2(posp[1], posp[0]), \
                    #      (xmax**2).sum()**0.5, math.atan2(xmax[1], xmax[0]) \

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
            ecc_max[i][j] =  max(eccplot)

    ecc_max_max = np.amax(ecc_max)
    ecc_max_min = np.amin(ecc_max)
    print "Maximum e of all runs:", ecc_max_max
    print "Minimum e of all runs:", ecc_max_min
    #plot colors/sizes
    color = np.where(ecc_max>0.5,'r','b')
    size = np.log(ecc_max/ecc_max_min) + 0.01

    #only need gto plot eccentricity
    fig = plt.figure()
    for i,mj in enumerate(Mj):
        for j,rj in enumerate(Rj):
            plt.scatter(mj, rj, s=size[i][j], c=color[i][j])
    plt.xlabel('$M_{j}$')
    plt.ylabel('$R_{j}$')
    plt.title('Outcomes')
    red_patch = mpatches.Patch(color='r', label='$ecc_{max} > 0.5$')
    blue_patch = mpatches.Patch(color='b', label='$ecc_{max} \leq 0.5$')
    plt.legend(handles=[red_patch, blue_patch], bbox_to_anchor=(1.04,1))

    #plt.tight_layout()
    #plt.show()
    fig.savefig("1b.png", bbox_inches='tight')

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
    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    #import pprint
    #pp = pprint.PrettyPrinter()
    #print 'Command-line parameters:'
    #pp.pprint(o.__dict__)
    main(**o.__dict__)
