import matplotlib.pyplot as plt
import numpy as np

#plt.style.use("science")

path = r"cmake-build-debug\Ablacio.txt"
matrix = []
z = np.linspace(0, 1, 100)
t = np.linspace(0, (0.49 * (0.01**2)) * 2000, 2000)
Z, T = np.meshgrid(z, t)

with open(path, "r") as f:
    lines = [line.rstrip() for line in f]

for line in lines:
    matrix.append(line.split(","))

for i in range(len(matrix)):
    for j in range(int(len(matrix[0]))):
        matrix[i][j] = float(matrix[i][j]) - 273.15

ax = plt.pcolormesh(Z, T, matrix)
cbar = plt.colorbar()
cbar.set_label("$T$ (ºC)")
plt.title("Temperatura")
plt.ylabel("$t$")
plt.xlabel("$z$")
plt.axvline(0.5 - 0.125, color="red")
plt.axvline(0.5 + 0.125, color="red")
plt.tight_layout()
plt.show()
