import odespy
import numpy as np
import matplotlib.pyplot as plt
import sys

def f(u, t):
    return -a*u

I = 1; a = 2; T = 6
dt = float(sys.argv[1]) if len(sys.argv) >= 2 else 0.75
Nt = int(round(T/dt))
t_mesh = np.linspace(0, Nt*dt, Nt+1)

solvers = [odespy.RK2(f),
           odespy.RK3(f),
           odespy.RK4(f),
           odespy.BackwardEuler(
               f, f_is_linear=True, jac=lambda u, t: -a)]

# (If f_is_linear is not specified in BackwardEuler, a
# nonlinear solver is invoked and this must be Newton
# (nonlinear_solver='Newton') since the default choice,
# nonlinear_solver='Picard', diverges in this problem.)

legends = []
for solver in solvers:
    solver.set_initial_condition(I)
    u, t = solver.solve(t_mesh)

    plt.plot(t, u)
    plt.hold('on')
    legends.append(solver.__class__.__name__)

# Compare with exact solution plotted on a very fine mesh
t_fine = np.linspace(0, T, 10001)
u_e = I*np.exp(-a*t_fine)
plt.plot(t_fine, u_e, '-') # avoid markers by specifying line type
legends.append('exact')

plt.legend(legends)
plt.title('Time step: %g' % dt)
plt.savefig('odespy1_dt_%g.png' % dt)
plt.savefig('odespy1_dt_%g.pdf' % dt)
plt.show()
