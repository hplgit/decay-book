import numpy as np

def differentiate(u, dt):
    dudt = np.zeros(len(u))
    for i in range(1, len(dudt)-1, 1):
        dudt[i] = (u[i+1] - u[i-1])/(2*dt)
    i = 0
    dudt[i] = (u[i+1] - u[i])/dt
    i = len(dudt)-1
    dudt[i] = (u[i] - u[i-1])/dt
    return dudt

def test_differentiate():
    """Test differentiate with a linear u."""
    # Expect exact results
    t = np.linspace(0, 4, 9)
    u = 2*t + 7
    dudt = differentiate(u, dt=t[1]-t[0])
    diff = abs(dudt - 2).max()
    tol = 1E-15
    assert diff < tol

def differentiate_vec(u, dt):
    dudt = np.zeros(len(u))
    dudt[1:-1] = (u[2:] - u[0:-2])/(2*dt)
    dudt[0] = (u[1] - u[0])/dt
    dudt[-1] = (u[-1] - u[-2])/dt
    #dudt[0] = (-3*u[0] + 4*u[1] - u[2])/(2*dt)
    #dudt[-1] = (3*u[-1] - 4*u[-2] + u[-3])/(2*dt)
    return dudt

def test_differentiate_vec():
    """Test differentiate_vec by comparing with differentiate."""
    t = np.linspace(0, 4, 9)
    u = 2*np.sin(t) + 7
    dudt_expected = differentiate(u, dt=t[1]-t[0])
    dudt_computed = differentiate_vec(u, dt=t[1]-t[0])
    diff = abs(dudt_expected - dudt_computed).max()
    tol = 1E-15
    assert diff < tol

def convergence(test_case):
    """Expect 2nd order, but what about the end points?"""
    # Use sympy to define and differentiate test examples
    import sympy as sym
    t = sym.symbols('t')
    # Symbolic expression and numerical interval
    if test_case == 1:
        u = sym.exp(-0.1*t)
        a = 0;  b = 10
    elif test_case == 2:
        u = sym.sqrt(t)
        a = 0.0001;  b = 2
    else:
        u = sym.sin(t)
        a = 0;  b = 4*np.pi

    dudt = sym.diff(u, t)
    print 'u:', u
    print 'dt/dt:', dudt, '\n'
    # Turn sympy expressions to numerical Python function
    u_func = sym.lambdify([t], u, modules='numpy')
    dudt_e = sym.lambdify([t], dudt, modules='numpy')

    print ' k      Nt   E/dt**2      E      Ei/dt**2    Ei'
    # Converge rate test, expect rate r=2, so chekc that errr/dt**r
    # is approx constant
    k_values = range(1, 20)
    for k in k_values:
        # Compute mesh
        Nt = 2**k;  dt = float(b-a)/Nt;  T = Nt*dt
        t = np.linspace(a, T, Nt+1)
        dudt = differentiate_vec(u_func(t), dt)
        # Approximate L2 integral of the error
        E  = np.sqrt(dt*np.sum((dudt_e(t) - dudt)**2))
        #E = np.abs(dudt_e(t).max() - dudt.max())
        Ei = np.sqrt(dt*np.sum((dudt_e(t)[1:-1] - dudt[1:-1])**2))
        #Ei = np.abs(dudt_e(t)[1:-1].max() - dudt[1:-1].max())
        print '%2d  %6d  %.2E  %.3E  %.2E  %.3E' % (k, Nt, E/dt**2, E, Ei/dt**2, Ei)
        #if k == 1:


#test_differentiate()
#test_differentiate_vec()
convergence(1)
convergence(2)
convergence(3)
