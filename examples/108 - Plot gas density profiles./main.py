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

    z           = np.linspace(-AU, AU)

    x           = rc / AU

    # disk_region
    plt.figure(figsize=(10, 11))

    plt.subplot(2, 2, 1)
    label = r"gas surface density $\Sigma_g(r)$  [kg m$^{-2}$]"
    plt.loglog(x, Sigma_g)
    plt.xlabel("distance from star $r$ [AU]")
    plt.ylabel(label)
    # plt.title(label)
    plt.grid()

    plt.subplot(2, 2, 2)
    label = r"gas midplane volume density $\rho_g^\text{mid}(r)$  [kg m$^{-3}$]"
    plt.loglog(x, rho_g_mid)
    plt.xlabel("distance from star $r$ [AU]")
    plt.ylabel(label)
    # plt.title(label)
    plt.grid()

    plt.subplot(2, 2, 3)
    label = r"gas midplane number density $N_g^\text{mid}(r)$  [m$^{-3}$]"
    plt.loglog(x, N_g_mid)
    plt.xlabel("distance from star $r$ [AU]")
    plt.ylabel(label)
    # plt.title(label)
    plt.grid()

    plt.subplot(2, 2, 4)
    label = r"gas volume density $N_g(r_0, z)$  [kg m$^{-3}$]"
    x = z / AU
    y = rho_g(r_0, z, rho_g_mid[i_0], H_p[i_0])
    plt.plot(x, y)
    plt.xlabel("height above/below midplane $z$ [AU]")
    plt.ylabel(label)
    # plt.title(label)
    plt.grid()

    plt.tight_layout()

    plt.savefig(Path(path_to_figures, "2x2 gas density.pdf"))
    plt.show()
    plt.close()

