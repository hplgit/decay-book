from numpy import *
from matplotlib.pyplot import *

def solver(I, a, b, T, dt, theta):
    """
    Solve u'=-a(t)*u + b(t), u(0)=I,
    for t in (0,T] with steps of dt.
    a and b are Python functions of t.
    """
    dt = float(dt)            # avoid integer division
    Nt = int(round(T/dt))     # no of time intervals
    T = Nt*dt                 # adjust T to fit time step dt
    u = zeros(Nt+1)           # array of u[n] values
    t = linspace(0, T, Nt+1)  # time mesh

    u[0] = I                  # assign initial condition
    for n in range(0, Nt):    # n=0,1,...,Nt-1
        u[n+1] = ((1 - dt*(1-theta)*a(t[n]))*u[n] + \
                  dt*(theta*b(t[n+1]) + (1-theta)*b(t[n])))/\
                  (1 + dt*theta*a(t[n+1]))
    return u, t

def test_linear_solution():
    """
    Test problem where u=c*t+I is the exact solution, to be
    reproduced (to machine precision) by any relevant method.
    """
    def u_exact(t):
        return c*t + I

    def a(t):
        return t**0.5  # can be arbitrary

    def b(t):
        return c + a(t)*u_exact(t)

    theta = 0.4; I = 0.1; dt = 0.1
    T = 4
    Nt = int(T/dt)  # no of steps

    c_values = [1E-5, 0.1, 1, 10, 100, 1000, 10000,
                1E+7, 1E+10, 1E+20]
    for c in c_values:
        u, t = solver(I=I, a=a, b=b, T=Nt*dt, dt=dt, theta=theta)
        u_e = u_exact(t)
        difference = abs(u_e - u).max()  # max deviation
        print 'c=%6g, difference=%g' % (c, difference)

test_linear_solution()
