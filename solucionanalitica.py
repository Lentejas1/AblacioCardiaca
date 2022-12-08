from numpy import sin, pi, exp

Tc = 273.15 + 36.5
k = 0.56
V = 40
cond = 0.472
P = cond * (V ** 2) / 2
N = 100


def f(x, t):
    suma = 0
    for n in range(1, 801):
        if n % 2 != 0:
            suma += 4 / (n * pi) * (1 - exp(-1 * t * (n * pi) ** 2)) / ((n * pi) ** 2) * sin(
                n * pi * x)
    return Tc * k / P + suma


lista_T = [round(f(i / N, 0.025) * (P / k) - 273.15, 5) for i in range(100)]
print(lista_T)
