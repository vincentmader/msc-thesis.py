import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from constants import AU
from models.solver import SolverV2
from models.disk import Disk, DiskRegion

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "108")
path_to_figures  = Path(PATH_TO_FIGURES, "108")
os.makedirs(path_to_figures, exist_ok=True)



def main():
    pass

def rho_g(r, z, rho_g_mid, H_p):
    return rho_g_mid * np.exp(- z**2 / 2 / H_p**2)


if __name__ == "__main__":
    cfg = Config()
    disk = Disk(cfg)
    disk_region = DiskRegion(cfg, disk)

    rc          = disk.rg.bin_centers
    Sigma_g     = disk.gas_surface_density
    rho_g_mid   = disk.midplane_gas_volume_density
    N_g_mid     = disk.midplane_gas_volume_number_density
    H_p         = disk.scale_height

    r_0         = disk_region.distance_to_star
    i_0         = disk.rg.index_from_value(r_0)

    z           = np.linspace(-1.5*AU, 1.5*AU)

    x           = rc / AU

    # disk_region
    plt.figure(figsize=(10, 11))

    plt.subplot(2, 2, 1)
    ylabel = r"Gas Surface Density $\Sigma_g(r)$  [kg m$^{-2}$]"
    label = r"Gas Surface Density $\Sigma_g(r)$"
    plt.loglog(x, Sigma_g, label=label)
    plt.scatter([r_0/AU], [Sigma_g[i_0]], marker="x", color="k", label="Distance from Star $r_0$", zorder=2)
    plt.xlabel("Distance from Star $r$ [AU]")
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.ylim(1e1, 1e5)
    # plt.title(label)
    plt.legend(loc="upper right")
    plt.grid()

    plt.subplot(2, 2, 2)
    ylabel = r"Midplane Gas Volume Density $\rho_g^\text{mid}(r)$  [kg m$^{-3}$]"
    label = r"Midplane Gas Volume Density $\rho_g^\text{mid}(r)$"
    plt.loglog(x, rho_g_mid, label=label)
    plt.scatter([r_0/AU], [rho_g_mid[i_0]], marker="x", color="k", label="Distance from Star $r_0$", zorder=2)
    plt.xlabel("Distance from Star $r$ [AU]")
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.ylim(1e-12, 1e-2)
    # plt.title(label)
    plt.legend(loc="upper right")
    plt.grid()

    plt.subplot(2, 2, 3)
    ylabel = r"Midplane Gas Number Density $N_g^\text{mid}(r)$  [m$^{-3}$]"
    label = r"Midplane Gas Number Density $N_g^\text{mid}(r)$"
    plt.loglog(x, N_g_mid, label=label)
    plt.scatter([r_0/AU], [N_g_mid[i_0]], marker="x", color="k", label="Distance from Star $r_0$", zorder=2)
    plt.xlabel("Distance from Star $r$ [AU]")
    plt.ylabel(ylabel)
    plt.xlim(x[0], x[-1])
    plt.ylim(1e15, 1e24)
    # plt.title(label)
    plt.legend(loc="upper right")
    plt.grid()

    plt.subplot(2, 2, 4)
    ylabel = r"Gas Volume Density $\rho_g(r_0, z)$  [kg m$^{-3}$]"
    label = r"Gas Volume Density $\rho_g(r_0, z)$"
    x = z / AU
    y = rho_g(r_0, z, rho_g_mid[i_0], H_p[i_0])
    plt.plot(x, y, label=label)
    plt.xlabel("Height above Midplane $z$ [AU]")
    plt.ylabel(label)
    plt.xlim(x[0], x[-1])
    plt.ylim(0, 8e-10)
    plt.legend(loc="upper right")
    # plt.title(label)
    plt.grid()

    plt.tight_layout()

    plt.savefig(Path(path_to_figures, "2x2 gas density.pdf"))
    plt.show()
    plt.close()

