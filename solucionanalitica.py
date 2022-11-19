import numpy as np


def b_n(m):
    return 4 / (m * np.pi)


def f(x, t):
    return sum([b_n(2 * n + 1) * (1 - np.exp(-2 * n + 1 ** 2 * np.pi ** 2 * t)) / (2 * n + 1 ** 2 * np.pi ** 2)
                * np.sin(2 * n + 1 * np.pi * x)] for n in range(100))


print(f(0.1, 1))
