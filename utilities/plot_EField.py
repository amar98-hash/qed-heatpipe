import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt(
    "data/Efield_map.dat"
)

row = data[:, 0].astype(int)
col = data[:, 1].astype(int)

Ex = data[:, 2]
Ey = data[:, 3]
Ez = data[:, 4]

nrows = row.max() + 1
ncols = col.max() + 1

Ex = Ex.reshape(nrows, ncols)
Ey = Ey.reshape(nrows, ncols)
Ez = Ez.reshape(nrows, ncols)

X, Y = np.meshgrid(
    np.arange(ncols),
    np.arange(nrows)
)

plt.quiver(
    X,
    Y,
    Ex,
    Ey
)

plt.axis("equal")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.show()


Emag = np.sqrt(
    Ex**2
    + Ey**2
    + Ez**2
)

plt.imshow(
    Emag,
    origin="lower"
)

plt.colorbar(label="|E| [V/m]")
plt.xlabel("x pixel")
plt.ylabel("y pixel")
plt.show()