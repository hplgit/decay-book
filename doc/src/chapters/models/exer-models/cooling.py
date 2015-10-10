import numpy as np

def cooling(T0, k, T_s, t_end, dt, theta=0.5):
    """
    Solve T'=-k(T-T_s(t)), T(0)=T0,
    for t in (0,t_end] with steps of dt.
    T_s(t) is a Python function of t.
    theta=0.5 means Crank-Nicolson, 1 is Backward
    Euler, and 0 is Forward Euler scheme.
    """
    dt = float(dt)                  # avoid integer division
    Nt = int(round(t_end/dt))       # no of time intervals
    t_end = Nt*dt                   # adjust to fit time step dt
    T = np.zeros(Nt+1)              # array of T[n] values
    t = np.linspace(0, t_end, Nt+1) # time mesh
    T[0] = T0                       # set initial condition
    for n in range(0, Nt):          # n=0,1,...,Nt-1
        T[n+1] = ((1 - dt*(1 - theta)*k)*T[n] + \
        dt*k*(theta*T_s(t[n+1]) + (1 - theta)*T_s(t[n])))/ \
        (1 + dt*theta*k)
    return T, t


def test_asymptotic():
    """
    Test that ``any'' initial condition leads to
    the same asymptotic behavior when T_s=constant.
    """
    import matplotlib.pyplot as plt
    plt.figure()
    T_s = 5.
    k = 1.2
    dt = 0.1
    tol = 0.01  # tolerance for testing asymptotic value
    t_end = 7    # make sure t_end is large enough for tol
    T0_values = [0, 2, 4, 5, 6, 8, 10] # test many cases

    for T0 in [0, 2, 4, 5, 6, 8, 10]:
        u, t = cooling(T0, k, lambda t: T_s, t_end, dt)
        plt.plot(t, u)

        assert abs(u[-1] - T_s) < tol, '%s != %s' % (u[-1], T_s)

    plt.legend(['T0=%g' % T0 for T0 in T0_values])
    plt.title('Testing asymptotic behavior T_s=%g' % T_s)
    plt.xlabel('t')
    plt.ylabel('T')
    plt.savefig('tmp1.png');  plt.savefig('tmp1.pdf')
    plt.show()


class Piecewise(object):
    """Class for holding a piecewise constant function."""
    def __init__(self, C0, C1, t_star):
        self.C0, self.C1 = C0, C1
        self.t_star = t_star

    def __call__(self, t):
        """
        Return value of piecewise constant function.
        t can be float or numpy array.
        """
        if isinstance(t, (float,int)):
            if t <= self.t_star:
                T_s = self.C0
            elif t > self.t_star:
                T_s = self.C1
        else:
            # assume numpy array
            T_s = np.piecewise(t,
                               [t <= self.t_star, t > self.t_star],
                               [self.C0, self.C1])
            # Alternative
            # T_s = np.where(t <= self.t_star, C0, C1)
        return T_s


def simulate_piecewise_constant_Ts():
    """
    Simulate scaled problem: T' = -(T - (2T_s-1)), T(0)=0,
    where T_s=1 for t < 3 and -0.5 for t > 3.
    """
    k = 1
    T0 = 0
    t_star = 3.0
    C0 = 1
    C1 = -0.5
    T_s = Piecewise(C0, C1, t_star)
    dt = t_star/100.0
    T, t = cooling(T0, k, T_s, t_end=3*t_star, dt=dt, theta=0.5)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(t, T)
    plt.xlabel('t');  plt.ylabel('u')
    plt.savefig('tmp2.png');  plt.savefig('tmp2.pdf')
    plt.show()

def T_exact_numint(self, t, k, T0, C0, C1, t_star):
    """
    Return exact (analytical) solution of the problem.
    The integral is the expression for the exact solution
    is approximated by an approximate trapezoidal rule
    (replacing weight 1/2 at the end points with just 1
    for simplicity).
    """
    T_s = Piecewise(C0, C1, t_star)
    t = np.linspace(0, t, 1000)  # integration points
    dt = t[1]-t[0]
    from numpy import exp  # for nicer math formula
    T = T0*exp(-k*t) + exp(-k*t)*k*dt*np.sum(
        self.T_s(t)*exp(k*t))
    return T

def evaluate_T_exact(t, k, T0, C0, C1, t_star, verbose=False):
    """
    Return exact (analytical) solution of the problem.
    Exact solution is produced by sympy.
    """
    exact0, exact1 = T_exact_symbolic()
    # exact0/1 works with t as numpy array
    if isinstance(t, (float,int)):
        if t < t_star:
            return exact0(t, C0, k, T0)
        else:
            return exact1(t, C0, C1, t_star, k, T0)
    else:
        # assume numpy array
        return np.where(
            t < t_star,
            exact0(t, C0, k, T0),
            exact1(t, C0, C1, t_star, k, T0))

def T_exact_symbolic(verbose=False):
    """Compute the exact solution formula via sympy."""
    # sol1: solution for t < t_star,
    # sol2: solution for t > t_star
    import sympy as sym
    T0 = sym.symbols('T0')
    k = sym.symbols('k', positive=True)
    # Piecewise linear T_sunction
    t, t_star, C0, C1 = sym.symbols('t t_star C0 C1')
    T_s = C0
    I = sym.integrate(sym.exp(k*t)*T_s, (t, 0, t))
    sol1 = T0*sym.exp(-k*t) + k*sym.exp(-k*t)*I
    sol1 = sym.simplify(sym.expand(sol1))
    if verbose:
        # Some debugging print
        print 'solution t < t_star:', sol1
        #print sym.latex(sol1)
    T_s = C1
    I = sym.integrate(sym.exp(k*t)*C0, (t, 0, t_star)) + \
        sym.integrate(sym.exp(k*t)*C1, (t, t_star, t))
    sol2 = T0*sym.exp(-k*t) + k*sym.exp(-k*t)*I
    sol2 = sym.simplify(sym.expand(sol2))
    if verbose:
        print 'solution t > t_star:', sol2
        #print sym.latex(sol2)

    # Convert to numerical functions
    exact0 = sym.lambdify([t, C0, k, T0],
                          sol1, modules='numpy')
    exact1 = sym.lambdify([t, C0, C1, t_star, k, T0],
                          sol2, modules='numpy')
    return exact0, exact1

def compare_numerical_and_exact_solution():
    """
    Compare exact and numerical solution with piecewise
    constant surrounding temperature. Use scaled problem
    from function simulate_piecewise_constant_Ts.
    """
    T0 = 0
    k = 1
    C0 = 1
    C1 = -0.5
    t_star = 3
    t_end = 7

    T_s = Piecewise(C0, C1, t_star)

    import matplotlib.pyplot as plt
    plt.figure()
    dt_values = [1, 0.5, 0.025]
    #dt_values = [0.025]
    for dt in dt_values:
        T, t = cooling(T0, k, T_s, t_end, dt, theta=0.5)
        plt.plot(t, T)

    t_e = np.linspace(0, t_end, 1001)  # find mesh for T_e
    # Could use sym.Rational(1,2) instead of 0.5, but not necessary
    # when we are not interested in symbolic formulas
    T_e = evaluate_T_exact(t_e, k, T0, C0, C1, t_star)
    plt.plot(t_e, T_e)
    plt.legend(['CN, dt=%g' % dt for dt in dt_values] + ['exact'])
    plt.title('T(t) for piecewise constant $T_s(t)$')
    plt.xlabel('t')
    plt.ylabel('T')
    plt.savefig('tmp3.png');  plt.savefig('tmp3.pdf')
    plt.show()


def test_discrete_solution():
    """
    Compare the numerical solution with an exact solution of the scheme
    when the T_s is constant.
    """
    T_s = 10
    T0 = 2
    k = 1.2
    dt = 0.1   # can use any mesh
    N_t = 6    # any no of steps will do
    t_end = dt*N_t
    t = np.linspace(0, t_end, N_t+1)

    for theta in [0, 0.5, 1, 0.2]:
        u, t = cooling(T0, k, lambda t: T_s , t_end, dt, theta)
        A = (1 - (1-theta)*k*dt)/(1 + theta*k*dt)
        u_discrete_exact = T_s + (T0-T_s)*A**(np.arange(len(t)))
        diff = np.abs(u - u_discrete_exact).max()
        print 'diff computed and exact discrete solution:', diff
        tol = 1E-14
        success = diff < tol
        assert success, 'diff=%g' % diff

if __name__=='__main__':
    #test_asymptotic()
    #simulate_piecewise_constant_Ts()
    #compare_numerical_and_exact_solution()
    test_discrete_solution()
