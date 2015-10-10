def estimate_k(T0, T1, Ts, dt):
    return float(T1 - T0)/(dt*(Ts - T0))

k = estimate_k(26.7, 25.8, 20, 3600)
print 'k=%g' % k

t = 0
dt = 1
if abs(k*dt) > 1:
    print("To large time step")

T = 37
Ts = 20
from cooling import cooling
while T > 25.8:
    T = T - k*dt*(T - Ts)
    t+= dt

minutes, seconds = divmod(t, 60)
hours, minutes = divmod(minutes, 60)
print """
The death occurred %d hours, %d minutes,
and %g seconds before 3am.""" % (hours, minutes, seconds)
