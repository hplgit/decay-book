import sympy as sym

# Define finite difference operators as functions

def D_f(u, dt, t):
    return (u(t + dt) - u(t))/dt

def D_b(u, dt, t):
    return (u(t) - u(t - dt))/dt

def D_c(u, dt, t):
    return (u(t + dt) - u(t - dt))/(2*dt)


def make_plot():
    def u(t):
        return sym.exp(-a*t)

    a, t, dt, p = sym.symbols('a t dt p')
    dudt = sym.diff(u(t), t)

    from numpy import logspace, exp
    from matplotlib.pyplot import (
        semilogx, legend, show, loglog, savefig)

    # Map operator function name to logical names
    operator2name = dict(
        D_f='forward', D_b='backward', D_c='central')
    legends = []
    for operator in D_f, D_b, D_c:
        E = operator(u, dt, t)/dudt
        # Expand, set p=a*dt, simplify
        E = sym.expand(E)
        E = E.subs(a*dt, p)
        E = sym.simplify(E)
        print '%s E:' % operator2name[operator.__name__], E
        print 'Taylor series:', E.series(p, 0, 3)
        latex_expr = sym.latex(E)

        E = sym.lambdify([p], E, modules='numpy')
        p_values = logspace(-6, -0.5, 101)
        y = E(p_values)
        semilogx(p_values, y)
        legends.append(operator2name[operator.__name__] +
                       ': $' + latex_expr + '$')
    legend(legends, loc='lower left')
    savefig('tmp.png'); savefig('tmp.pdf')
    show()

make_plot()
