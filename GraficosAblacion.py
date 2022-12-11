import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, exp, sin

plt.style.use("science")

Tc = 273.15 + 36.5
k = 0.56
V = 40
cond = 0.472
P = cond * (V ** 2) / 2
dz = 0.01
Z = 2
z = np.linspace(0, dz * 100 * Z, 101)
alpha = 0.56 / (3683 * 1081)


def f(x, t):
    suma = 0
    for n in range(1, 801):
        if n % 2 != 0:
            suma += 4 / (n * pi) * (1 - exp(-1 * t * (n * pi) ** 2)) / ((n * pi) ** 2) * sin(n * pi * x)
    return Tc * k / P + suma


def plot_resultado(path_f, M_f, N_f, ratio, show=False):
    matrix = []
    dt = ratio * dz ** 2
    T = 0.02 ** 2 / alpha
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
    if show:
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
    plt.savefig("figures/resultat_grafic.png", dpi=300)
    if show:
        plt.show()
    print(f"t_max={careful * dt * T} s")
    show = False


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

def tiempo_optimo():
    T = 0.02 ** 2 / alpha
    T_f = (36.5 + 273.15) * k / P
    T_opt = (50 + 273.15) * k / P
    t = 0
    while T_f < T_opt:
        t += 1E-6
        T_f = f(0.375, t)

    return t*T

def plot_error(ratio, metodo, path_f, filename, output_data=False, show=False):
    array_numerico = t0025file(path_f)
    if metodo == "explicito":
        dt = ratio * dz ** 2
    if metodo == "implicito":
        dt = ratio * dz
    M = 0.025 // dt + 1
    lista_E_T = [abs(float(array_numerico[i]) - f(i / 100, (M - 1) * dt) * (P / k) + 273.15) for i in range(101)]
    plt.figure(figsize=(4, 3))
    plt.scatter(z, lista_E_T, color="royalblue", s=1)
    plt.xlabel("$z$ (cm)")
    plt.ylabel("$E_{T}$ ($^\circ$C)")
    mean = np.mean(lista_E_T)
    std = np.std(lista_E_T)
    plt.axhline(mean, linestyle="--", color="r", label=r"$\bar{E_T}=$" + "{:.2g} $^\circ$C".format(mean))
    plt.tight_layout()
    plt.legend()
    if output_data:
        print(metodo, " ratio =", ratio, " M =", M, "\n")
        print(f"Media error = {np.mean(lista_E_T)} K")
        print("std =", std, " K")
        print("-" * 20)
    filename += ".png"
    plt.savefig(filename, dpi=300)
    if show:
        plt.show()
    show = False


def plot_analitico(show=False):
    lista_T = [f(i / 100, 0.025) * (P / k) - 273.15 for i in range(101)]
    plt.figure(figsize=(4, 3))
    plt.plot(z, lista_T, color="royalblue")
    plt.xlabel("$z$ (cm)")
    plt.ylabel("$T$ ($^\circ$C)")
    plt.tight_layout()
    plt.axvline(0.75, color="red")
    plt.axvline(1.25, color="red")
    plt.savefig("figures/Analítico_ta.png", dpi=300)
    if show:
        plt.show()
    show = False


def ta0025_plot(file, label, figurefilename, show=False, **kwargs):
    plt.figure(figsize=(5, 3))
    plt.scatter(z, t0025file(file),
                label=label,
                s=0.1)
    plt.xlabel("$z$ (m)")
    plt.ylabel("$T$ ($^\circ$C)")
    plt.legend(loc='lower center')
    plt.tight_layout()
    try:
        plt.scatter(z, t0025file(kwargs["file2"]), label=(kwargs["label2"]), s=0.1)
    except KeyError:
        pass

    plt.legend(loc='lower center')
    plt.tight_layout()
    plt.savefig(figurefilename, dpi=300)
    if show:
        plt.show()
    show = False


def tramesa():
    with open("tramesa/tramesa_explicit.txt", "w") as file:
        file.write("Explícito 0.25")
        file.write(str(t0025file("data/Ablacio_Explicit_025.txt")))
        array = t0025file("data/Ablacio_Explicit_025.txt")
        file.write(str([abs(float(array[i]) - f(i / 100, (0.025 // (0.25 * dz ** 2) * 0.25 * dz ** 2) * (P / k) +
                                                273.15)) for i in range(101)]))
        file.write("Explícito 0.49")
        file.write(str(t0025file("data/Ablacio_Explicit_049.txt")))
        array = t0025file("data/Ablacio_Explicit_049.txt")
        file.write(str([abs(float(array[i]) - f(i / 100, (0.025 // (0.49 * dz ** 2) * 0.49 * dz ** 2) * (P / k) +
                                                273.15)) for i in range(101)]))
        file.write("Explícito 0.51")
        file.write(str(t0025file("data/Ablacio_Explicit_051.txt")))
        array = t0025file("data/Ablacio_Explicit_051.txt")
        file.write(str([abs(float(array[i]) - f(i / 100, (0.025 // (0.51 * dz ** 2) * 0.51 * dz ** 2) * (P / k) +
                                                273.15)) for i in range(101)]))

    with open("tramesa/tramesa_implicit.txt", "w") as file:
        file.write("Implícito 0.50")
        file.write(str(t0025file("data/Ablacio_Implicit050.txt")))
        array = t0025file("data/Ablacio_Implicit050.txt")
        file.write(str([abs(float(array[i]) - f(i / 100, (0.025 // (0.5 * dz ** 2) * 0.5 * dz ** 2) * (P / k) +
                                                273.15)) for i in range(101)]))
        file.write("Implícito 1.00")
        file.write(str(t0025file("data/Ablacio_Implicit100.txt")))
        array = t0025file("data/Ablacio_Implicit100.txt")
        file.write(str([abs(float(array[i]) - f(i / 100, (0.025 // (1 * dz ** 2) * 1 * dz ** 2) * (P / k) +
                                                273.15)) for i in range(101)]))


# GRAFICA EL ANALÍTICO
print(f"Tiempo óptimo = {tiempo_optimo()}")
plot_analitico()

# GRÁFICOS A t_a=0.025 (por defecto solo los guarda, show=True para enseñarlo)
ta0025_plot(file="data/Ablacio_Explicit_025.txt", file2="data/Ablacio_Explicit_049.txt",
            label="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.25$",
            label2="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.49$",
            figurefilename="figures/t_aexp025049.png", show=False)
ta0025_plot(file="data/Ablacio_Explicit_051.txt", label="$\dfrac{\Delta \hat{t}}{(\Delta \hat{z})^2}=0.51$",
            figurefilename="figures/t_aexp051.png", show=False)
ta0025_plot(file="data/Ablacio_Implicit050.txt", file2="data/Ablacio_Implicit100.txt",
            label="$\dfrac{\Delta \hat{t}}{\Delta \hat{z}}=0.50$", label2="$\dfrac{\Delta \hat{t}}{\Delta \hat{z}}=1$",
            figurefilename="figures/t_aimp050100.png", show=False)

# GRÁFICOS ERRORES A t_0.025
plot_error(0.25, metodo="explicito", path_f="data/Ablacio_Explicit_025.txt", filename="figures/Error_T_exp_0.25_ta",
           output_data=True, show=False)
plot_error(0.49, metodo="explicito", path_f="data/Ablacio_Explicit_049.txt", filename="figures/Error_T_exp_0.49_ta",
           output_data=True, show=False)
plot_error(0.51, metodo="explicito", path_f="data/Ablacio_Explicit_051.txt", filename="figures/Error_T_exp_0.51_ta",
           output_data=True, show=False)
plot_error(0.50, metodo="implicito", path_f="data/Ablacio_Implicit050.txt", filename="figures/Error_T_imp_0.50_ta",
           output_data=True, show=False)
plot_error(1, metodo="implicito", path_f="data/Ablacio_Implicit100.txt", filename="figures/Error_T_imp_1.00_ta",
           output_data=True, show=False)

# GRAFICA EL MEJOR MÉTODO Y SU EVOLUCIÓN TEMPORAL Y EL t_optimo
plot_resultado("data/Ablacio_Explicit_ResultatFinal.txt", 2000, 101, 0.25, show=False)

# CREA LOS ARCHIVOS REQUERIDOS EN EL APARTADO DE ENVÍO DEL CV
tramesa()

