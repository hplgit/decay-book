# Import from ../../genz-src
import sys, os
sys.path.insert(
    0, os.path.join(os.pardir, os.pardir, 'genz', 'src-genz'))
from decay_vc import solver

# Scaled model: u' = -u + b(t), b=1 for t < 1 else 0
eps, t = solver(I=0, a=lambda t: 1,
                b=lambda t: 1 if t < 1 else 0,
                T=4, dt=0.01, theta=0.5)

import matplotlib.pyplot as plt
plt.plot(t, eps)
plt.xlabel('Dimensionless time'); plt.ylabel('10*strain')
plt.savefig('tmp.png');  plt.savefig('tmp.pdf')
plt.show()
