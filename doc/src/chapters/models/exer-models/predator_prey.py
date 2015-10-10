import odespy
import numpy as np
#import scitools.std as plt
import matplotlib.pyplot as plt
fig_counter = 0  # used in plot file names

def simulate(alpha, beta, gamma, delta, T=20):
    def f(u, t, alpha, beta, gamma):
        H, L = u
        return [alpha*H - L*H, beta*L*H - gamma*L]


    t = np.linspace(0, T, 1501)
    solver = odespy.RK4(f, f_args=(alpha, beta, gamma))
    solver.set_initial_condition([1, delta])
    u, t = solver.solve(t)
    H = u[:,0]
    L = u[:,1]

    plt.figure()
    global fig_counter
    fig_counter += 1
    plt.plot(t, H, t, L)
    plt.legend(['H', 'L'])
    plt.title(r'$\alpha=%g$, $\beta=%g$, $\gamma=%g$, $\delta=%g$'
              % (alpha, beta, gamma, delta))
    plt.savefig('tmp%d.png' % fig_counter)
    plt.savefig('tmp%d.pdf' % fig_counter)
    return H, L, t

def demo():
    simulate(alpha=1.5, beta=0.5, gamma=0.8, delta=0.5, T=20)
    simulate(alpha=0.5, beta=0.5, gamma=0.2, delta=0.5, T=50)

def simulate2(mu, nu, omega, T=20):
    def f(u, t, mu):
        H, L = u
        return [H*(1-L), mu*L*(H-1)]


    t = np.linspace(0, T, 1501)
    solver = odespy.RK4(f, f_args=(mu,))
    solver.set_initial_condition([nu, omega])
    u, t = solver.solve(t)
    H = u[:,0]
    L = u[:,1]

    plt.figure()
    global fig_counter
    fig_counter += 1
    plt.plot(t, H, t, L)
    plt.legend(['H', 'L'])
    plt.title(r'$\mu=%g$, $\nu=%g$, $\omega=%g$' % (mu, nu, omega))
    plt.savefig('tmp%d.png' % fig_counter)
    plt.savefig('tmp%d.pdf' % fig_counter)
    return H, L, t

def demo2():
    simulate2(mu=0.5, nu=0.5, omega=1.5, T=20)
    simulate2(mu=2, nu=0.5, omega=1.5, T=20)

if __name__ == '__main__':
    #demo()
    demo2()
    plt.show()
    raw_input()
