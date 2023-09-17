import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import KernelAxis
    from config import Config
    from dust import particle_mass_from_radius
    from kernel import Kernel
    from kernel.mass_conservation import test_mass_conservation
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
    from visualization.kernel.mass_conservation import KernelMassConservationSubplot
except ModuleNotFoundError as e:
    raise e


cfg = Config()
kernel = Kernel(cfg)

mg = kernel.mg
mc = mg.grid_cell_centers
ac = mg.particle_radii

rho_s = cfg.dust_particle_density

def plot_1():
    s1 = PcolorMatrixSubplot(
        ac, ac, kernel.K_gain, 
        title="kernel gain contribution $G_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        symmetrized=True,
        z_limits=(1e-20, 1e-7),
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, -kernel.K_loss,
        title="kernel loss contribution $L_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        symmetrized=True,
        z_limits=(1e-20, 1e-7),
    )
    p = GridspecPlot([s1, s2], add_slider=True)
    p.render()


def plot_2():
    # TODO Move the below elsewhere? Best: In `KernelMassConservationSubplot` definition.
    sum_ij = test_mass_conservation(mg, kernel.K)
    def custom_format_coord(x, y):
        x, y = particle_mass_from_radius(x, rho_s), particle_mass_from_radius(y, rho_s)
        x, y = mg.index_from_value(x), mg.index_from_value(y)
        i, j = int(x), int(y)  # TODO Index Convention? (irrelevant due to symmetry)
        text = f"sum_k K_kij = {sum_ij[i, j]:.2}, {i = }, {j = }"
        return text

    s = KernelMassConservationSubplot(
        mg, kernel.K,
        title=r"kernel mass error $\sum_k m_k K_{kij}$",
        axis=KernelAxis.Radius,
    )
    p = GridspecPlot([s])
    p.axes[1].format_coord = custom_format_coord
    p.render()


def main():
    plot_1()
    plot_2()
