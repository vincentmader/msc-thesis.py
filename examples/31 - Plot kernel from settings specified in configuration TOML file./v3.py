import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_radius_from_mass, particle_mass_from_radius
    from kernel import Kernel
    from kernel.mass_conservation import test_mass_conservation
    from visualization.kernel.v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.v3_2023_08_14.gridspec_plot import GridspecPlot
    from visualization.kernel.v4_2023_09_13.mass_conservation import KernelMassConservationSubplot
except ModuleNotFoundError as e:
    raise e


cfg = Config()
kernel = Kernel(cfg)

mg = kernel.mg
mc = mg.grid_cell_centers
rho_s = cfg.dust_particle_density
ac = particle_radius_from_mass(mc, rho_s)

def plot_1():
    s1 = PcolorMatrixSubplot(
        ac, ac, kernel.K_gain, 
        title="kernel gain contribution $G_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        symmetrize=True,
        z_limits=(1e-20, 1e-7),
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, -kernel.K_loss,
        title="kernel loss contribution $L_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        symmetrize=True,
        z_limits=(1e-20, 1e-7),
    )
    p = GridspecPlot([s1, s2], add_slider=True)
    p.render()


def plot_2():
    sum_ij = test_mass_conservation(cfg, mg, kernel.K)

    def custom_format_coord(x, y):
        x, y = particle_mass_from_radius(x, rho_s), particle_mass_from_radius(y, rho_s)
        x, y = mg.index_from_value(x), mg.index_from_value(y)
        i, j = int(x), int(y)  # TODO Index Convention? (irrelevant due to symmetry)
        text = f"sum_k K_kij = {sum_ij[i, j]:.2}, {i = }, {j = }"
        return text

    s = KernelMassConservationSubplot(
        kernel,
        title=r"kernel mass error $\sum_k m_k K_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
    )
    p = GridspecPlot([s])
    p.axes[1].format_coord = custom_format_coord
    p.render()


def main():
    plot_1()
    plot_2()

# def format_coord(self, x, y):
#     return self.custom_format_coord(x, y)
