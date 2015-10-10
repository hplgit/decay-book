from cooling import cooling
from numpy import pi, sin

def T_s(t):
    return Tm + a*sin((2*pi/P)*t)

Tm = 25
a = 2.5
P_values = [3600, 600, 3600*6]
k = 0.05/60
T0 = 5

import matplotlib.pyplot as plt
for P in P_values:
    T, t = cooling(T0, k, T_s, t_end=8*3600, dt=600/40)
    plt.plot(t, T)
plt.plot(t, T_s(t), 'k--')  # T_s for largest P to show amplitude
legends = ['P=1 h', 'P=10 min', 'P=6 h', '$T_s$']
plt.legend(legends, loc='lower right')
plt.xlabel('t'); plt.ylabel('T')
plt.savefig('tmp.png');  plt.savefig('tmp.pdf')
plt.show()
