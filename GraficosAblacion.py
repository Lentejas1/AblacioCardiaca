import matplotlib.pyplot as plt
import numpy as np

plt.style.use("science")


def plot_explicit(path_f, M_f, N_f, ratio):
    matrix = []
    dz = 0.01
    dt = ratio * dz ** 2
    alpha = 0.56 / (3683 * 1081)
    Z = 2
    T = 0.02 ** 2 / alpha
    z = np.linspace(0, dz * 100 * Z, 100)
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
    plt.savefig("implicito5heatmap.png")
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
    plt.savefig("implicito5grafico.png")
    plt.show()

    print(f"t_max={careful * dt * T} s")


M = 2000
N = 100


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
    plt.savefig("implicito1heatmap.png")
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
    plt.savefig("implicito1grafico.png")
    plt.show()

    print(f"t_max={careful * dt * T} s")


N = 100

plot_implicit(r"cmake-build-debug\Ablacio_Implicit1.txt", 6, N, 1)
