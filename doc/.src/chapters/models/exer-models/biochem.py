import odespy
import numpy as np
import matplotlib.pyplot as plt
import sys

def solver(alpha, beta, epsilon, T, dt=0.1):
    def f(u, t):
        Q, P, S, E = u
        return [
            alpha*(E*S - Q),
            beta*Q,
            -E*S + (1-beta/alpha)*Q,
            (-E*S + Q)/epsilon,
            ]

    Nt = int(round(T/dt))
    t_mesh = np.linspace(0, Nt*dt, Nt+1)

    solver = odespy.RK4(f)
    solver.set_initial_condition([0, 0, 1, 1])
    u, t = solver.solve(t_mesh)
    Q = u[:,0]
    P = u[:,1]
    S = u[:,2]
    E = u[:,3]
    return Q, P, S, E

def demo():
    alpha = 1.5
    beta = 1
    epsilon = 0.1
    T = 8
    dt = float(sys.argv[1]) if len(sys.argv) >= 2 else 0.01
    Q, P, S, E = solver(alpha, beta, epsilon, T, dt)
    plt.plot(t, Q, t, P, t, S, t, E)
    plt.legend(['complex', 'product', 'substrate', 'enzyme'],
               loc='upper right')
    plt.title('alpha=%g, beta=%g, epsilon=%g' %
              (alpha, beta, epsilon))
    plt.savefig('tmp.png');  plt.savefig('tmp.pdf')
    plt.show()

if __name__ == '__main__':
    demo()
