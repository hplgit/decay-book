from dc_leapfrog import *

def leapfrog(I, a, b, T, dt, filter=True,
             first_step='ForwardEuler',
             gamma=0.6):
    dt = float(dt)
    N = int(round(T/dt))  # no of intervals
    u = np.zeros(N+1)
    t = np.linspace(0, N*dt, N+1)

    u[0] = I
    if first_step == 'ForwardEuler':
        u[1] = (1 - dt*a(t[0])*u[0] + dt*b(t[0])
    elif first_step == 'CrankNicolson':
        u[1] = ((1 - 0.5*a(t[0])*dt)*u[0] + \
                0.5*dt*(b(t[0]) + b(t[1])))/\
               (1 + 0.5*dt*a(t[1]))
    for n in range(1, N):
        u[n+1] = u[n-1] - 2*dt*a(t[n])*u[n] + 2*dt*b(t[n])
        if filtered:
            u[n] = u[n] + gamma*(u[n-1] - 2*u[n] + u[n+1])
    return u, t

def main():
    pass

if __name__ == '__main__':
    main()
