import numpy as np

def mesh_function(f, t):
    u = np.zeros(len(t))  # or t.size
    for i in range(len(t)):
        u[i] = f(t[i])
    return u

def demo():
    def f(t):
        if t <= 3:
            return np.exp(-t)
        else:
            return np.exp(-3*t)

    # Compute mesh and mesh function
    t = np.linspace(0, 4, 41)
    u = mesh_function(f, t)

    # Plot
    import matplotlib.pyplot as plt
    plt.plot(t, u)
    plt.xlabel('t')
    plt.ylabel('mesh function')
    plt.savefig('tmp.png'); plt.savefig('tmp.pdf')
    plt.show()
