from functions1a import *


def main(N, seed, eps, dt, t_end, v0, eta):
    if eps <= 0.0: eps = tiny
    
    # Initial conditions.

    #cluster 1
    mass1,pos1,vel1 = initialize(N, seed, v0)
    pos1 = pos1 + np.array([2.0, 0.5, 0])
    vel1 = vel1 + np.array([-0.5, 0.0,0])

    #cluster 2
    mass2,pos2,vel2 = initialize(N, seed, v0)
    pos2 = pos2 + np.array([-2.0, -0.5, 0])
    vel2 = vel2 + np.array([0.5, 0.0,0])

    #combinte them
    mass = np.concatenate((mass1, mass2), axis=0)
    pos = np.concatenate((pos1, pos2), axis=0)
    vel = np.concatenate((vel1, vel2), axis=0)

    #move to center of mass frame
    m = mass.reshape(N*2,1)/mass.sum()
    pos -= (m*pos).sum(axis=0)
    vel -= (m*vel).sum(axis=0)

    t = 0.0
    dt0 = dt
    timestep = 0
    acc,tau = acceleration2(mass, pos, vel, eps**2)
    dt = eta*dt0*tau

    time_interval = np.array([1.0*x for x in range(int(math.ceil(t_end+3)))])
    ind = 1

    # Initial diagnostics.
    E0 = energy(mass, pos, vel, eps**2)
    print 'Initial E =', E0
    output(t, E0, mass, pos, vel, eps**2, timestep)

    # Run forward to specified time.
    tplot = []
    dEplot = []

    Eplot = []
    Uplot = []
    Rplot = []
    KEplot = []
    
    while t < t_end-0.5*dt:
        #apply variable time step with eta
        t,pos,vel,tau = step(t, mass, pos, vel, eps**2, dt)
        timestep += 1
        dt = eta*dt0*tau#next dt

        E = energy(mass, pos, vel, eps**2)
        
        #a,e,Erel,h = orbital_elements(mass[0], mass[1], pos[0],
        #                              pos[1], vel[0], vel[1], eps**2)

        U = potential_energy(mass, pos, eps**2)
        R = virial_rad(mass, pos, eps**2)
        KE = kinetic_energy(mass, vel)
        tplot.append(t)
        dEplot.append(E-E0)

        Eplot.append(E)
        Uplot.append(U)
        Rplot.append(R)
        KEplot.append(KE)

        #print energy each time unit
        if t>time_interval[ind]-dt and t<=time_interval[ind]:
            print "t = ", t
            print "dE = ",abs(E-E0)
            ind +=1

        if (abs(E-E0) > 1e-4):
            print "Energy error too great, breaking loop"
            print "dE = ", E-E0
            break



    # Final diagnostics.

    output(t, E0, mass, pos, vel, eps**2, timestep)

    print "timesteps taken: ", timestep

    #plotting
    plt.plot(tplot, Eplot, label="Total Energy")
    plt.plot(tplot, Uplot, label="Potential Energy")
    plt.plot(tplot, KEplot, label="Kinetic Energy")
    plt.plot(tplot, Rplot, label="Virial Radius")
    plt.xlabel('t')
    plt.ylabel('E')
    plt.title('$\eta = %s$' % eta)
    plt.legend(bbox_to_anchor=(1.04,1))
    plt.savefig("1a.png", bbox_inches='tight')
    plt.show()
    plt.clf()


def new_option_parser():
    from optparse import OptionParser
    result = OptionParser()
    result.add_option("-n",
                      dest="N", type="int", default ="50",
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
                      dest="t_end", type="float", default ="30.0",#"30.0",
                      help="integration interval [%default]")
    result.add_option("-v",
                      dest="v0", type="float", default ="0.25",
                      help="initial 2-body v [%default]")
    result.add_option("-b",
                      dest="eta", type="float", default ="10.0",
                      help="max energy error allowed [%default]")

    return result

if __name__ in ('__main__'):
    o, arguments  = new_option_parser().parse_args()
    import pprint
    pp = pprint.PrettyPrinter()
    print 'Command-line parameters:'
    pp.pprint(o.__dict__)
    main(**o.__dict__)
