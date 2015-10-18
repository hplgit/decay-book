from numpy import *
from matplotlib.pyplot import *

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

def explore(I, a, T, dt_values, theta=0.5):
    """
    Run cases with the solver, compute error measure,
    and plot the error.
    """
    figure() # create new plot
    for dt in dt_values:
        print 'dt', dt
        u, t = solver(I, a, T, dt, theta)
        # Numerical solution
        u_e = u_exact(t, I, a)
        e = u_e - u
        plot(t, e)
        hold('on')
    legend(['dt=%g' % dt for dt in dt_values])
    xlabel('t')
    ylabel('u')
    title('theta=%g' % theta)
    theta2name = {0: 'FE', 1: 'BE', 0.5: 'CN'}
    savefig('%s_%g.png' % (theta2name[theta], dt))
    savefig('%s_%g.pdf' % (theta2name[theta], dt))

def main(I, a, T, dt_values, theta_values=(0, 0.5, 1)):
    for theta in theta_values:
        explore(I, a, T, dt_values, theta)

dt = 0.1
dt_values = [dt, dt/4, dt/8]
main(I=1, a=2, T=5, dt_values=dt_values)
show()
