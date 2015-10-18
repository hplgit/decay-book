import sympy as sym

T0 = sym.symbols('T0')
k = sym.symbols('k', positive=True)
# Piecewise linear function
t, t_a, C0, C1 = sym.symbols('t t_a C0 C1')
T_s = C0
I = sym.integrate(sym.exp(k*t)*T_s, (t, 0, t))
sol1 = T0*sym.exp(-k*t) + k*sym.exp(-k*t)*I
sol1 = sym.simplify(sym.expand(sol1))
print sol1
print sym.latex(sol1)
T_s = C1
I = sym.integrate(sym.exp(k*t)*C0, (t, 0, t_a)) + \
    sym.integrate(sym.exp(k*t)*C1, (t, t_a, t))
sol2 = T0*sym.exp(-k*t) + k*sym.exp(-k*t)*I
sol2 = sym.simplify(sym.expand(sol2))
print sol2
print sym.latex(sol2)

# Try a sine variation
T_s = sym.sin(2*t)
I = sym.integrate(sym.exp(k*t)*T_s, (t, 0, t))
sol3 = T0*sym.exp(-k*t) + k*sym.exp(-k*t)*I
print sol3
sol3 = sym.lambdify([t, k, T0], sol3, modules='numpy')

# Convert to numerical functions

sol1 = sym.lambdify([t, C0, k, T0], sol1, modules='numpy')
sol2 = sym.lambdify([t, C0, C1, t_a, k, T0], sol2, modules='numpy')

def T_exact(t, C0, C1, t_a, k, T0):
    if isinstance(t, (float,int)):
        return sol1(t, C0, k, T0) if t < t_a \
               else sol2(t, C0, C1, t_a, k, T0)
    else:
        # Vectorized version
        return np.where(t < t_a, sol1(t, C0, k, T0),
                        sol2(t, C0, C1, t_a, k, T0))

import matplotlib.pyplot as plt
import numpy as np
t = np.linspace(0, 8, 801)
#T = T_exact(t=t, C0=1, C1=2, t_a=3, k=1, T0=0)
T = T_exact(t=t, C0=1, C1=-2, t_a=3, k=1, T0=0)
plt.plot(t, T)
plt.title('Sudden change in surrounding temperature')

plt.figure()
plt.plot(t, sol3(t=t, k=0.5, T0=0), t, sol3(t=t, k=2, T0=0), t, np.sin(2*t))
plt.legend(['$T(t;k=0.5)$', '$T(t;k=2)$', '$T_s(t)=\sin 2t$'])
plt.title('Sinusoidal surrounding temperature')
plt.show()
