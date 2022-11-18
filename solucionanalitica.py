import numpy as np


def b_n(m):
    return 4 / (m * np.pi)


def f(x, t):
    s = 0
    for n in range(100):
        if n != 0:
            m = 2 * n + 1
            s += b_n(m) * (1 - np.exp(-m ** 2 * np.pi ** 2 * t)) / (m ** 2 * np.pi ** 2) * np.sin(m * np.pi * x)
    return s


print(f(0.1, 1))
