from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
import numpy as np

N_phi = 100
N_r   = 100
R_min = 1e-2
R_max = 1e-1  # Maximum radius of the disk

radial_coords    = np.linspace(R_min, R_max, num=N_r)
azimuthal_coords = np.linspace(0, 2 * np.pi, num=N_phi)

r, theta = np.meshgrid(radial_coords, azimuthal_coords)
surface_density = 1.0 / r

x = r * np.cos(theta)
y = r * np.sin(theta)

# ═════════════════════════════════════════════════════════════════════════════

plt.figure(figsize=(8, 8))

contour = plt.contourf(x, y, surface_density, cmap="Blues", levels=500)
plt.colorbar(contour, label="Gas Surface Density", norm=LogNorm())

plt.title("Protoplanetary Disk - Gas Surface Density")
plt.axis("equal")
plt.xlim(-R_max, R_max)
plt.ylim(-R_max, R_max)
plt.gca().spines["right"].set_visible(False)
plt.gca().spines["top"].set_visible(False)
plt.gca().spines["bottom"].set_visible(False)
plt.gca().spines["left"].set_visible(False)
plt.xticks([])
plt.yticks([])

plt.show()
