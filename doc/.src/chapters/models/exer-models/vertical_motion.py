# HPL: To reuse this file as a module, all computations should be
# encapuslated in function (i.e., the time stepping loop and the plotting
# in your case). Otherwise this is a very clear, clean and to-the-point
# solution. One of the very best I've seen. Despite the simple math,
# programming this exercise has resulted in a lot of messy codes...

import sys
import matplotlib.pyplot as plt
from math import pi,log10
from numpy import zeros, linspace

try:
	N = int(sys.argv[1])
	T = float(sys.argv[2])
except:
	print'Please provide commandline arguments as follows:'
	sys.exit()

#------ initialization-----------
g = 9.81 # m/s^2
d= 0.5*10**(-3) #diameter in meters
dt = T/N
mu = 1.8*10**(-5) # viscosity
rho = 0.0#1.2 # density (of air)
rho_b = 1000 #mass density of body
C_d = 0.45 #drag coefficient (for a sphere)
V = (4./3)*pi*(d/2)**3 # volume of body [m**3]
A = pi*(d/2)**2	#cross sectional area of body perpendicular to motion [m**2]
a2 = (C_d*A*rho)/(2*rho_b*V)
a1 = (3*pi*d*mu)/(rho_b*V)
b = g*(rho/rho_b -1)
Re = 0
v = zeros(N)

#------- functions ------------------
def Reynolds (v,d,rho,mu):
	return (rho*d*v)/mu
def constants(Re):
	return a1 if Re<1.0 else a2
def step(Re,v,a):
	step_stokes = ((1+0.5*a*dt)*v +b*dt)/(1+0.5*dt*a)
	step_quadratic = (v+dt*b)/(1+dt*a*abs(v))
	return step_stokes if Re<1.0 else step_quadratic

#-------- calculations ---------------
for i in range(N-1):
	Re = Reynolds(abs(v[i]),d,rho,mu)
	v[i+1] = step(Re,v[i],constants(Re))

#-------- verification ----------------
v_ver = [-g*i for i in linspace(0,T,N)]
error = [(v[i]-v_ver[i])/v_ver[i] for i in range(len(v))]
print 'for the straight line approx, the error is in the %.6f th term'% log10(max(error))
#print error
plt.plot(v)
plt.show()

