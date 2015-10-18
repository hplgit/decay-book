import decay
import nose.tools as nt
import numpy as np

def u_discrete_exact(n, I, a, theta, dt):
    """Return exact discrete solution of the numerical schemes."""
    dt = float(dt)  # avoid integer division
    A = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)
    return I*A**n

def test_u_discrete_exact():
    """
    Compare result from solver against
    formula for the discrete solution.
    """
    theta = 0.8; a = 2; I = 0.1; dt = 0.8
    Nt = int(8/dt)  # no of steps
    u, t = decay.solver(I=I, a=a, T=Nt*dt, dt=dt, theta=theta)
    u_de = np.array(
        [u_discrete_exact(n, I, a, theta, dt)
         for n in range(Nt+1)])
    diff = np.abs(u_de - u).max()
    nt.assert_almost_equal(diff, 0, delta=1E-14)
    tol = 1E-14
    success = diff < tol
    assert success

def test_potential_integer_division():
    """Choose variables that can trigger integer division."""
    theta = 1; a = 1; I = 1; dt = 2
    Nt = 4
    u, t = decay.solver(I=I, a=a, T=Nt*dt, dt=dt, theta=theta)
    u_de = np.array(
        [u_discrete_exact(n, I, a, theta, dt)
         for n in range(Nt+1)])
    diff = np.abs(u_de - u).max()
    tol = 1E-14
    nt.assert_almost_equal(diff, 0, delta=1E-14)

def test_read_command_line_positional():
    # Decide on a data set of input parameters
    I = 1.6;  a = 1.8;  T = 2.2;  theta = 0.5
    dt_values = [0.1, 0.2, 0.05]
    # Expected return from read_command_line_positional
    expected = [I, a, T, theta, dt_values]
    # Construct corresponding sys.argv array
    sys.argv = [sys.argv[0], str(I), str(a), str(T), 'CN'] + \
               [str(dt) for dt in dt_values]
    computed = decay.read_command_line_positional()
    for expected_arg, computed_arg in zip(expected, computed):
        assert expected_arg == computed_arg

def test_read_command_line_argparse():
    I = 1.6;  a = 1.8;  T = 2.2;  theta = 0.5
    dt_values = [0.1, 0.2, 0.05]
    # Expected return from read_command_line_positional
    expected = [I, a, T, theta, dt_values]
    # Construct corresponding sys.argv array
    cml = '%s --a %s --I %s --T %s --scheme CN --dt ' % \
          (sys.argv[0], a, I, T)
    cml = cml + ' '.join([str(dt) for dt in dt_values])
    sys.argv = cml.split()
    computed = decay.read_command_line_argparse()
    for expected_arg, computed_arg in zip(expected, computed):
        assert expected_arg == computed_arg
