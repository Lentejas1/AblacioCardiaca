import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, exp, sin

plt.style.use("science")

Tc = 273.15 + 36.5
k = 0.56
V = 40
cond = 0.472
P = cond * (V ** 2) / 2


def f(x, t):
    suma = 0
    for n in range(1, 801):
        if n % 2 != 0:
            suma += 4 / (n * pi) * (1 - exp(-1 * t * (n * pi) ** 2)) / ((n * pi) ** 2) * sin(n * pi * x)
    return Tc * k / P + suma


def plot_explicit(path_f, M_f, N_f, ratio):
    matrix = []
    dz = 0.01
    dt = ratio * dz ** 2
    alpha = 0.56 / (3683 * 1081)
    Z = 2
    T = 0.02 ** 2 / alpha
    z = np.linspace(0, dz * 100 * Z, 101)
    t = np.linspace(0, dt * M_f * T, M_f)
    with open(path_f, "r") as f:
        lines = [line.rstrip() for line in f]

    for line in lines:
        matrix.append(line.split(","))

    for i in range(len(matrix)):
        for j in range(int(len(matrix[0]))):
            matrix[i][j] = float(matrix[i][j]) - 273.15

    plt.figure(figsize=(4, 3))
    plt.pcolormesh(z, t, matrix)
    cbar = plt.colorbar()
    cbar.set_label("$T$ ($^\circ$C)")
    plt.ylabel("$t$ (s)")
    plt.xlabel("$z$ (cm)")
    plt.axvline(0.75, color="red")
    plt.axvline(1.25, color="red")
    plt.tight_layout()
    plt.savefig("figures/resultat_heatmap.png", dpi=300)

    # CHEQUEO MAX Y MUERTE
    unhealthy_cured = []
    healthy_unsafe = []
    careful_bool = False
    alive = True
    for i in range(len(matrix)):
        unhealthy_counter_cured = 0
        healthy_counter_unsafe = 0
        for j in range(int(len(matrix[0]))):
            if (3 / 8 * N_f <= j) and (j <= 5 / 8 * N_f):
                if float(matrix[i][j]) >= 50:
                    if float(matrix[i][j]) <= 80:
                        unhealthy_counter_cured += 1
                    elif alive:
                        alive = False
                        dead = i
            else:
                if float(matrix[i][j]) > 50:
                    if not careful_bool:
                        careful_bool = True
                        careful = i
                    healthy_counter_unsafe += 1
        unhealthy_cured.append(unhealthy_counter_cured)
        healthy_unsafe.append(healthy_counter_unsafe)

    plt.figure(figsize=(4, 3))

    plt.plot(t, unhealthy_cured, label="Células curadas")
    plt.plot(t, healthy_unsafe, label="Células sanas muertas")
    plt.xlabel("$t$ (s)")
    plt.ylabel("$n$")
    if careful_bool:
        plt.axvline(careful * dt * T, color="red")
    if not alive:
        plt.axvline(dead, label="Trombosis", color="red")
    plt.legend(loc='upper left')
    plt.xlim(careful * dt * T - 20, careful * dt * T + 20)
    plt.tight_layout()
    plt.savefig("figures/resultat_grafic.png", dpi=300)

    print(f"t_max={careful * dt * T} s")


def plot_implicit(path_f, M_f, N_f, ratio):
    matrix = []
    dz = 0.01
    dt = ratio * dz
    alpha = 0.56 / (3683 * 1081)
    Z = 2
    T = 0.02 ** 2 / alpha
    z = np.linspace(0, dz * 100 * Z, 100)
    t = np.linspace(0, dt * M_f * T, M_f)
    print(t)

    with open(path_f, "r") as f:
        lines = [line.rstrip() for line in f]

    for line in lines:
        matrix.append(line.split(","))

    for i in range(len(matrix)):
        for j in range(int(len(matrix[0]))):
            matrix[i][j] = float(matrix[i][j]) - 273.15

    plt.figure(figsize=(4, 3))
    plt.pcolormesh(z, t, matrix)
    cbar = plt.colorbar()
    cbar.set_label("$T$ ($^\circ$C)")
    plt.ylabel("$t$ (s)")
    plt.xlabel("$z$ (cm)")
    plt.axvline(0.75, color="red")
    plt.axvline(1.25, color="red")
    plt.tight_layout()
    plt.savefig("figures/implicito1heatmap.png")
    plt.show()

    # CHEQUEO MAX Y MUERTE
    unhealthy_cured = []
    healthy_unsafe = []
    careful_bool = False
    alive = True
    for i in range(len(matrix)):
        unhealthy_counter_cured = 0
        healthy_counter_unsafe = 0
        for j in range(int(len(matrix[0]))):
            if (3 / 8 * N_f <= j) and (j <= 5 / 8 * N_f):
                if float(matrix[i][j]) >= 50:
                    if float(matrix[i][j]) <= 80:
                        unhealthy_counter_cured += 1
                    elif alive:
                        alive = False
                        dead = i
            else:
                if float(matrix[i][j]) > 50:
                    if not careful_bool:
                        careful_bool = True
                        careful = i
                    healthy_counter_unsafe += 1
        unhealthy_cured.append(unhealthy_counter_cured)
        healthy_unsafe.append(healthy_counter_unsafe)

    plt.figure(figsize=(4, 3))

    plt.plot(t, unhealthy_cured, label="Células curadas")
    plt.plot(t, healthy_unsafe, label="Células sanas muertas")
    plt.xlabel("$t$ (s)")
    plt.ylabel("$n$")
    if careful_bool:
        plt.axvline(careful * dt * T, color="red")
    if not alive:
        plt.axvline(dead, label="Trombosis", color="red")
    plt.legend(loc='upper left')
    plt.xlim(careful * dt * T - 20, careful * dt * T + 20)
    plt.tight_layout()
    plt.savefig("figures/implicito1grafico.png")
    plt.show()

    print(f"t_max={careful * dt * T} s")


def t0025file(path_f):
    matrixt00025 = []
    arrayimp00025 = []
    with open(path_f, "r") as f:
        lines = [line.rstrip() for line in f]

    for line in lines:
        matrixt00025.append(line.split(","))

    for i in range(len(matrixt00025)):
        for j in range(int(len(matrixt00025[0]))):
            matrixt00025[i][j] = float(matrixt00025[i][j]) - 273.15

    for i in range(101):
        arrayimp00025.append(matrixt00025[-1][i])

    return arrayimp00025


def plot_error(path_f, filename):
    array_numerico = t0025file(path_f)
    lista_E_T = [abs(float(array_numerico[i]) - f(i / 100, 0.025) * (P / k) + 273.15) for i in range(101)]
    plt.figure(figsize=(4, 3))
    plt.scatter(z, lista_E_T, color="royalblue", s=1)
    plt.xlabel("$z$ (cm)")
    plt.ylabel("$E_{T}$ ($^\circ$C)")
    mean = np.mean(lista_E_T)
    std = np.std(lista_E_T)
    plt.axhline(mean, linestyle="--", color="r", label=r"$\bar{E_T}=$" + "{:.2g} $^\circ$C".format(mean))
    plt.tight_layout()
    plt.legend()
    # print(f"Media={np.mean(lista_E_T)} K")
    # print("std=", std)
    filename += ".png"
    plt.savefig(filename, dpi=300)
    # plt.show()


def plot_analitico():
    lista_T = [f(i / 100, 0.025) * (P / k) - 273.15 for i in range(101)]
    plt.figure(figsize=(4, 3))
    plt.scatter(z, lista_T, color="royalblue", s=1)
    plt.xlabel("$z$ (cm)")
    plt.ylabel("$T$ ($^\circ$C)")
    plt.tight_layout()
    plt.axvline(0.75, color="red")
    plt.axvline(1.25, color="red")
    plt.savefig("figures/Analítico_ta.png", dpi=300)


dz = 0.01
Z = 2
z = np.linspace(0, dz * 100 * Z, 101)

plt.figure(figsize=(5, 3))
plt.scatter(z, t0025file("cmake-build-debug/Ablacio_Explicit025.txt")
            , label="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.25$", s=0.1)
plt.scatter(z, t0025file("cmake-build-debug/Ablacio_Explicit049.txt")
            , label="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.49$", s=0.1)
plt.xlabel("$z$ (m)")
plt.ylabel("$T$ ($^\circ$C)")
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("figures/t_aexp025049.png", dpi=300)

plt.figure(figsize=(5, 3))
plt.scatter(z, t0025file("cmake-build-debug/Ablacio_Explicit051.txt"),
            label="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.51$",
            s=0.1)
plt.xlabel("$z$ (m)")
plt.ylabel("$T$ ($^\circ$C)")
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("figures/t_aexp050100.png", dpi=300)

plt.figure(figsize=(5, 3))
plt.scatter(z, t0025file("cmake-build-debug/Ablacio_Implicit050.txt"),
            label="$\dfrac{\Delta \hat{t}}{\Delta \hat{z}}=0.5$", s=0.1)
plt.scatter(z, t0025file("cmake-build-debug/Ablacio_Implicit100.txt"),
            label="$\dfrac{\Delta \hat{t}}{\Delta \hat{z}}=1$", s=0.1)
plt.xlabel("$z$ (m)")
plt.ylabel("$T$ ($^\circ$C)")
plt.legend(loc='lower center')
plt.tight_layout()
plt.savefig("figures/t_aimp050100.png", dpi=300)

plot_explicit("cmake-build-debug/Ablacio_Explicit_ResultatFinal.txt", 2000, 101, 0.25)

plot_error("cmake-build-debug/Ablacio_Explicit025.txt", "figures/Error_T_exp_0.25_ta")
plot_error("cmake-build-debug/Ablacio_Explicit049.txt", "figures/Error_T_exp_0.49_ta")
plot_error("cmake-build-debug/Ablacio_Explicit051.txt", "figures/Error_T_exp_0.51_ta")
plot_error("cmake-build-debug/Ablacio_Implicit050.txt", "figures/Error_T_imp_0.50_ta")
plot_error("cmake-build-debug/Ablacio_Implicit100.txt", "figures/Error_T_imp_1.00_ta")
plot_analitico()

