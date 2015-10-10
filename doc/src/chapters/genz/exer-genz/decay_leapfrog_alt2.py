# Here is a set of more reusable functions doing the same
# steps as above, but more modularized

# Global symbols to be seen by all local functions
t, I, dt, a, b, c = sym.symbols('t I dt a b c')

def ode_source_term(u):
    """Return the terms in the ODE that the source term
    must balance, here u' + a*u
    u is symbolic Python function of t."""
    return sym.diff(u(t), t) + a*u(t)

def residual_discrete_eq(u):
    """Return the residual of the discrete eq. with u inserted."""
    R = D2t(u, dt) + a*u(t) - b
    return sym.simplify(R)

def residual_discrete_eq_step1(u):
    """Return the residual of the discrete eq. at the first
    step with u inserted."""
    R = Dtp(u, dt) + a*u(t) - b  # Forward Euler
    R = R.subs(t, 0)  # t=0 in the rhs of the first step eq.
    return sym.simplify(R)

def D2t(u, dt):
    """Return 2nd-order finite difference for u_t.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - u(t-dt))/(2*dt)

def Dtp(u, dt):
    """Return 1st-order forward finite difference for u_t.
    u is a symbolic Python function of t.
    """
    return (u(t+dt) - u(t))/dt

def main(u):
    """
    Given some chosen solution u (as a function of t, implemented
    as a Python function), use the method of manufactured solutions
    to compute the source term f, and check if u also solves
    the discrete equations.
    """
    print '=== Testing exact solution: %s ===' % u(t)
    print "Initial condition u(0)=%s" % u(t).subs(t, 0)

    # Method of manufactured solution requires fitting b

    global b
    b = ode_source_term(u)
    b = sym.simplify(b)

    # Residual in discrete equations (should be 0)
    print 'residual step1:', residual_discrete_eq_step1(u)
    print 'residual:', residual_discrete_eq(u)

def linear():
    def u_e(t):
        """Return chosen linear exact solution."""
        # General linear function u_e = c*t + d
        # Initial condition u(0)=I requires d=I
        return c*t + I

    main(u_e)

def quadratic():
    main(lambda t: t**2 + c*t + I)

#linear
quadratic()
