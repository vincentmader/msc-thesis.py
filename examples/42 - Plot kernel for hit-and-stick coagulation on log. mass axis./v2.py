import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_radius_from_mass
    from kernel import Kernel
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
except ModuleNotFoundError as e:
    raise e


cfg = Config(
    enable_fragmentation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)
kernel_1 = Kernel(cfg)

mg = kernel_1.mg
mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

s1 = PcolorMatrixSubplot(
    ac, ac, kernel_1.K_gain, 
    title="kernel gain contribution $G_{kij}$",
    xlabel="particle radius $a_j$ [m]",
    ylabel="particle radius $a_i$ [m]",
    symmetrized=True,
)
s2 = PcolorMatrixSubplot(
    ac, ac, -kernel_1.K_loss,
    title="kernel loss contribution $L_{kij}$",
    xlabel="particle radius $a_j$ [m]",
    symmetrized=True,
)

def main():
    p = GridspecPlot([s1, s2], add_slider=True)
    p.render()
