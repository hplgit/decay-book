import numpy as np

# Note: this implementation does not have b(t) and a(t)

def backward2step(I, a, T, dt, first_step='BackwardEuler'):
    """2nd-order backward scheme."""
    dt = float(dt)
    N = int(round(T/dt))  # no of intervals
    u = np.zeros(N+1)
    t = np.linspace(0, N*dt, N+1)

    u[0] = I
    if first_step == 'CrankNicolson':
        u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]  # Crank-Nicolson 1. step
    elif first_step == 'BackwardEuler':
        u[1] = 1/(1 + dt*a)*u[0]                   # Backward Euler 1. step
    for n in range(1, N):
        u[n+1] = (4*u[n] - u[n-1])/(3 + 2*dt*a)
    return u, t

