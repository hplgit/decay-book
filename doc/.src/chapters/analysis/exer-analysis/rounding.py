def forward(u, t, dt):
    return (u(t+dt) - u(t))/dt

def backward(u, t, dt):
    return (u(t) - u(t-dt))/dt

def centered(u, t, dt):
    return (u(t+dt) - u(t-dt))/(2*dt)

def differentiation():
    from math import exp

    def u(t):
        return exp(-t)

    def dudt(t):
        return -u(t)

    print ' k    forward   backward  centered'
    t = 0.5
    for k in range(0, 63, 4):
        dt = 2**(-k)
        error_forward  = abs(dudt(t) - forward (u, t, dt))
        error_backward = abs(dudt(t) - backward(u, t, dt))
        error_centered = abs(dudt(t) - centered(u, t, dt))
        print '%2d    %.2E  %.2E  %.2E' % \
              (k, error_forward, error_backward, error_centered)

def trapezoidal(f, a, b, n):
    """
    Trapezoidal rule for integral of f from a to b
    using n intervals.
    """
    h = float(b - a)/n
    I = 0.5*(f(a) + f(b))
    for i in range(1, n):
        I += f(a + i*h)
    return h*I

def midpoint(f, a, b, n):
    """
    Midpoint-rule for integral of f from a to b
    using n intervals.
    """
    h = float(b - a)/n
    I = 0
    for i in range(0, n):
        I += f(a + (i+0.5)*h)
    return h*I

def integration():
    from math import exp

    def u(t):
        return exp(-t)

    def U(a, b):
        """Integral of u(t) from a to b."""
        return (-u(b) - (-u(a)))

    print ' k    trapezoidal    midpoint'
    a = 0
    b = 4
    for k in range(1, 23, 2):
        n = 2**k
        error_trapezoidal  = abs(U(a, b) - trapezoidal(u, a, b, n))
        error_midpoint = abs(U(a, b) - midpoint(u, a, b, n))
        print '%2d     %.2E      %.2E' % \
              (k, error_trapezoidal, error_midpoint)

differentiation()
integration()
