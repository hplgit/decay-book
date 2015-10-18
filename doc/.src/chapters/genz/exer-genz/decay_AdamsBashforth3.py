import numpy as np

# Note: this implementation does not have b(t) and a(t)

def Adams_Bashforth_3(I, a, T, dt):
    """
    Third-order Adams-Bashforth scheme (4 time levels), with
    Crank-Nicolson 2nd order scheme for the first two steps.
    """
    dt = float(dt)
    N = int(round(T/dt))  # no of intervals
    u = np.zeros(N+1)
    t = np.linspace(0, N*dt, N+1)

    u[0] = I
    u[1] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[0]  # Crank-Nicolson 1. step
    u[2] = (1 - 0.5*a*dt)/(1 + 0.5*dt*a)*u[1]  # Crank-Nicolson 2. step
    for n in range(2, N):
        u[n+1] = u[n] - dt/12.*a*(23*u[n] - 16*u[n-1] + 5*u[n-2])
    return u, t
