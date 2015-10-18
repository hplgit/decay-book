import numpy as np
import matplotlib.pyplot as plt

def solver(I, a, u_crit, dt, theta):
    """
    Solve u'=-a*u, u(0)=I, for t in (0,t_m] until u <= u_crit
    with steps of dt. Return t_m.
    """
    # Use list for u and t since we do not know how many points
    # that are needed
    dt = float(dt)               # avoid integer division
    u = []
    t = []

    u.append(I)                  # assign initial condition
    t.append(0)
    while u[-1] > u_crit:
        u_new = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[-1]
        u.append(u_new)
        t.append(t[-1] + dt)
    return t[-1]


half_life = 5730
a  = np.log(2)/half_life
print 'Age:', solver(I=1, a=a, u_crit=0.084, dt=10, theta=0.5)

half_life_min = 5730 - 40
half_life_max = 5730 + 40
a_min = np.log(2)/half_life_min
a_max = np.log(2)/half_life_max
age_min = solver(I=1, a=a_max, u_crit=0.084, dt=10, theta=0.5)
age_max = solver(I=1, a=a_min, u_crit=0.084, dt=10, theta=0.5)
print 'Uncertainty: [%g, %g]' % (age_min, age_max)

