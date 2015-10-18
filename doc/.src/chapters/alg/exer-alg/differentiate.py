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

test_differentiate()
test_differentiate_vec()
