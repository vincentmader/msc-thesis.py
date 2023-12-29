import matplotlib.pyplot as plt
import numpy as np
import random
from tqdm import tqdm


N_phi = 100
N_r   = 100
R_min = 1e-2
R_max = 1e+5  # Maximum radius of the disk

radial_coords    = np.linspace(R_min, R_max, num=N_r)
azimuthal_coords = np.linspace(0, 2 * np.pi, num=N_phi)

r, theta = np.meshgrid(radial_coords, azimuthal_coords)
surface_density = 1.0 / r

# ═════════════════════════════════════════════════════════════════════════════

x = []
y = []

p = 1 / radial_coords
for r, n_r in tqdm(zip(radial_coords, p)):
    for i in range(int(n_r)):
        X = 100
        phi = random.randint(0, int(2*np.pi * X)) / X

        x.append(r * np.cos(theta))
        y.append(r * np.sin(theta))
        plt.scatter([x], [y], label="")

print(len(x))
plt.show()
