from numpy import *

# Exercise a

def solver(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    dt = float(dt)            # avoid integer division
    Nt = int(round(T/dt))     # no of time intervals
    T = Nt*dt                 # adjust T to fit time step dt
    u = zeros(Nt+1)           # array of u[n] values
    t = linspace(0, T, Nt+1)  # time mesh

    u[0] = I                  # assign initial condition
    for n in range(0, Nt):    # n=0,1,...,Nt-1
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t

def u_exact(t, I, a):
    return I*exp(-a*t)

def numerical_and_exact(theta, I, a, T, dt):
    """Compare the numerical and exact solution in a plot."""
    u, t = solver(I=I, a=a, T=T, dt=dt, theta=theta)

    t_e = linspace(0, T, 1001)        # fine mesh for u_e
    u_e = u_exact(t_e, I, a)
    return u, t, u_e, t_e

def demo(dt):
    from matplotlib.pyplot import (
        plot, xlabel, ylabel, legend, title, savefig, show)
    for theta in [0, 0.5, 1]:
        u, t, u_e, t_e = numerical_and_exact(
            I=1, a=-1, T=2.5, dt=dt, theta=theta)
        xlabel('t')
        ylabel('u')
        plot(t,   u)

    plot(t_e, u_e, 'k-')  # black line
    legend(['FE', 'CN', 'BE', 'exact'], loc='upper left')
    title('Timestep: %g' % dt)
    savefig('tmp_%g.png' % dt); savefig('tmp_%g.pdf' % dt)
    show()

# Exercise b

def plot_amplification_factors(names):
    # Substitute -p by p since a is negative for a growth model

    def A_exact(p):
        return exp(p)

    def A(p, theta):
        return (1+(1-theta)*p)/(1-theta*p)

    def amplification_factor(names):
        # Use SciTools since it adds markers to colored lines
        from scitools.std import (
            plot, title, xlabel, ylabel, hold, savefig,
            axis, legend, grid, show, figure)
        figure()
        curves = {}
        p = linspace(0, 3, 99)
        curves['exact'] = A_exact(p)
        plot(p, curves['exact'])
        hold('on')
        name2theta = dict(FE=0, BE=1, CN=0.5)
        for name in names:
            curves[name] = A(p, name2theta[name])
            plot(p, curves[name])
            axis([p[0], p[-1], -20, 20])
            #semilogy(p, curves[name])
        plot([p[0], p[-1]], [0, 0], '--')  # A=0 line
        title('Amplification factors')
        grid('on')
        legend(['exact'] + names, loc='lower left', fancybox=True)
        xlabel(r'$p=-a\cdot dt$')
        ylabel('Amplification factor')
        savefig('A_growth.png'); savefig('A_growth.pdf')
        #show()

    amplification_factor(names)

if __name__ == '__main__':
    import sys
    dt_values = [float(arg) for arg in sys.argv[1:]]
    for dt in dt_values:
        demo(dt)  # Exercise a
    plot_amplification_factors('FE BE CN'.split()) # Exercise b
