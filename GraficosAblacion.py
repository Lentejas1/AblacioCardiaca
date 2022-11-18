import matplotlib.pyplot as plt
import numpy as np

#plt.style.use("science")

M = 1000

path = r"cmake-build-debug\Ablacio.txt"
matrix = []
z = np.linspace(0, 1, 100)
t = np.linspace(0, (0.49 * (0.01**2)) * M, M)

with open(path, "r") as f:
    lines = [line.rstrip() for line in f]
    print(len(lines))

for line in lines:
    matrix.append(line.split(","))
print(len(matrix))

for i in range(len(matrix)):
    for j in range(int(len(matrix[0]))):
        matrix[i][j] = float(matrix[i][j]) - 273.15

print(len(matrix))
ax = plt.pcolormesh(z, t, matrix)
cbar = plt.colorbar()
cbar.set_label("$T$ (ºC)")
plt.title("Temperatura")
plt.ylabel("$t$")
plt.xlabel("$z$")
plt.axvline(0.5 - 0.125, color="red")
plt.axvline(0.5 + 0.125, color="red")
plt.tight_layout()
plt.show()
