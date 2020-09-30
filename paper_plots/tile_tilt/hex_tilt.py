import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.linalg import lstsq

# Elevation, East, North
hex = np.array(
    [
        ["HexS01", 376.542, -9.985, 158.452],
        ["HexS02", 376.411, 4.006, 158.456],
        ["HexS03", 376.270, 18.022, 158.439],
        ["HexS04", 376.158, 32.013, 158.442],
        ["HexS05", 376.596, -17.013, 146.309],
        ["HexS06", 376.428, -2.985, 146.316],
        ["HexS07", 376.323, 11.022, 146.311],
        ["HexS08", 376.160, 25.002, 146.320],
        ["HexS09", 376.042, 39.014, 146.304],
        ["HexS10", 376.578, -23.966, 134.201],
        ["HexS11", 376.416, -9.983, 134.185],
        ["HexS12", 376.314, 4.005, 134.206],
        ["HexS13", 376.164, 17.990, 134.182],
        ["HexS14", 376.056, 32.009, 134.184],
        ["HexS15", 375.935, 46.016, 134.221],
        ["HexS16", 376.602, -30.991, 122.075],
        ["HexS17", 376.457, -16.989, 122.070],
        ["HexS18", 376.332, -2.972, 122.058],
        ["HexS19", 376.089, 25.015, 122.055],
        ["HexS20", 375.961, 39.014, 122.090],
        ["HexS21", 375.822, 53.019, 122.052],
        ["HexS22", 376.510, -23.994, 109.937],
        ["HexS23", 376.345, -9.980, 109.944],
        ["HexS24", 376.229, 3.996, 109.951],
        ["HexS25", 376.088, 18.025, 109.959],
        ["HexS26", 375.978, 32.024, 109.944],
        ["HexS27", 375.834, 46.022, 109.949],
        ["HexS28", 376.369, -16.987, 97.810],
        ["HexS29", 376.228, -2.998, 97.819],
        ["HexS30", 376.149, 11.019, 97.833],
        ["HexS31", 375.997, 25.004, 97.836],
        ["HexS32", 375.917, 39.021, 97.820],
        ["HexS33", 376.251, -9.987, 85.702],
        ["HexS34", 376.169, 4.018, 85.698],
        ["HexS35", 376.049, 18.018, 85.716],
        ["HexS36", 375.920, 32.013, 85.713],
    ]
)


Z = hex[:, 1].astype(float)
E = hex[:, 2].astype(float)
N = hex[:, 3].astype(float)

Z = Z - min(Z)
E = E - min(E)
N = N - min(N)


tmp_A = []
tmp_b = []
for i in range(len(E)):
    tmp_A.append([E[i], N[i], 1])
    tmp_b.append(Z[i])
b = np.matrix(tmp_b).T
A = np.matrix(tmp_A)
fit, residual, rnk, s = lstsq(A, b)

X, Y = np.meshgrid(np.linspace(min(E), max(E), 256), np.linspace(min(N), max(N), 256))
ele = np.zeros(X.shape)
for r in range(X.shape[0]):
    for c in range(X.shape[1]):
        ele[r, c] = fit[0] * X[r, c] + fit[1] * Y[r, c] + fit[2]

# Move origin to bottom left instead of top left
ele = np.flip(ele, axis=0)
X = np.flip(X, axis=0)
Y = np.flip(Y, axis=0)

i_min = np.unravel_index(np.argmin(ele, axis=None), ele.shape)
i_max = np.unravel_index(np.argmax(ele, axis=None), ele.shape)

fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(111)
ax.set_aspect("equal")
im = ax.imshow(
    ele, extent=(min(E) - 2, max(E) + 2, min(N) - 2, max(N) + 2), cmap="viridis",
)
plt.scatter(
    E, N, c=Z, marker="s", cmap="viridis", s=77, edgecolor="black", linewidth=0.4
)
plt.xlabel("East [m]")
plt.ylabel("North [m]")
plt.title(
    f"MinMax Coords (m): [({X[i_min]:.3f}E, {Y[i_min[::-1]]:.3f}N, {ele[i_min]:.3f}Z), ({X[i_max]:.3f}E, {Y[i_max[::-1]]:.3f}N, {ele[i_max]:.3f}Z)]",
    fontsize=11,
)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.1)
fig.colorbar(im, cax=cax)
plt.tight_layout()
#  plt.show()
plt.savefig("./interp_maps/beam_offsets/hex_tilt.png")
