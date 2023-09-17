import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, KernelAxis
    from config import Config
    from dust import particle_radius_from_mass, particle_mass_from_radius
    from kees_kernel import create_coag_kernel
    from kernel import Kernel
    from kernel.mass_conservation import test_mass_conservation
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
    from visualization.kernel.kernel import KernelSubplot
    from visualization.kernel.mass_conservation import KernelMassConservationSubplot
    from visualization.v1.mass_conservation import KernelMassConservationPlot
except ModuleNotFoundError as e:
    raise e


def plot_1(
    scale: str,
    mg: DiscreteMassAxis,
    Ks: tuple[np.ndarray, np.ndarray, np.ndarray],
):
    K_vinc, K_kees, K_compare = Ks
    axis = KernelAxis.Bin if scale == "lin" else KernelAxis.Radius

    s1 = KernelSubplot(
        mg, K_vinc, 
        title="$K_{kij}^{vinc}$",
        axis=axis,
        symmetrized=True,
        scales=(scale, scale, "lin"),
        z_limits=(-1, 1),
        cmap="bwr",
    )
    s2 = KernelSubplot(
        mg, K_kees, 
        title="$K_{kij}^{kees}$",
        axis=axis,
        symmetrized=True,
        scales=(scale, scale, "lin"),
        z_limits=(-1, 1),
        cmap="bwr",
    )
    s3 = KernelSubplot(
        mg, K_compare, 
        title=r"$\Delta K_{kij}=K_{kij}^{kees}-K_{kij}^{vinc}$",
        axis=axis,
        symmetrized=True,
        scales=(scale, scale, "lin"),
        z_limits=(-1, 1),
        cmap="bwr",
    )
    p = GridspecPlot([s1, s2, s3], add_slider=True)
    p.render()

def plot_2(
    scale: str,
):
    # if scale == "lin":
    #     xlabel = "particle radius $a_j$ [m]",
    #     ylabel = "particle radius $a_i$ [m]",
    # else:
    #     xlabel = "particle radius $a_j$ [m]",
    #     ylabel = "particle radius $a_i$ [m]",

    # axis = KernelAxis.Bin if scale == "lin" else KernelAxis.Radius
    # s = KernelSubplot(mg, kernel_vinc.K, axis=axis)
    # p = GridspecPlot([s], add_slider=True)
    # p.render()

    z = test_mass_conservation(cfg, mg, K_vinc)
    s1 = PcolorMatrixSubplot(
        x, y, z,
        title=r"kernel mass error $\sum_k m_k K_{kij}$",
        # xlabel="particle radius $a_j$ [m]",
        # ylabel="particle radius $a_i$ [m]",
        scales=(scale, scale, "log"),
    )
    z = test_mass_conservation(cfg, mg, K_kees)
    s2 = PcolorMatrixSubplot(
        x, y, z,
        title=r"kernel mass error $\sum_k m_k K_{kij}$",
        # xlabel="particle radius $a_j$ [m]",
        # ylabel="particle radius $a_i$ [m]",
        scales=(scale, scale, "log"),
    )
    p = GridspecPlot([s1, s2])

    # sum_ij = test_mass_conservation(cfg, mg, kernel_vinc.K)
    # def custom_format_coord(x, y):
    #     x, y = particle_mass_from_radius(x, rho_s), particle_mass_from_radius(y, rho_s)
    #     x, y = mg.index_from_value(x), mg.index_from_value(y)
    #     i, j = int(x), int(y)  # TODO Index Convention? (irrelevant due to symmetry)
    #     text = f"sum_k K_kij = {sum_ij[i, j]:.2}, {i = }, {j = }"
    #     return text

    # p.axes[1].format_coord = custom_format_coord
    p.render()

def plot_3(cfg, mg, K):
    KernelMassConservationPlot(cfg, mg, K).show()


setups = [
    ("lin", (2   , 52  , 50)),
    ("log", (1e-4, 1e+4, 50)),
]
if __name__ == "__main__":
    for scale, (m_min, m_max, N_m) in setups:
        cfg = Config(
            mass_axis_scale=scale, mass_resolution=N_m,
            mass_min_value=m_min, mass_max_value=m_max,
            enable_coagulation=True, enable_fragmentation=False,
            enable_physical_collisions=False, relative_velocity_components=[],
            enable_collision_sampling=False, # If true: Supply `ijs`!
            enable_cancellation_handling=False,
        )
    
        kernel_vinc = Kernel(cfg)
        K_vinc, mg = kernel_vinc.K, kernel_vinc.mg
        mc, ac = mg.grid_cell_centers, mg.particle_radii
    
        R_coll = np.ones(shape=[N_m] * 2) # TODO Redefine `R` (and `K` ?)
        K_kees = create_coag_kernel(mc, R_coll)
        
        K_vinc_sym = np.array([0.5 * (z + z.T) for z in K_vinc])
        K_diff = K_kees - K_vinc_sym
        # K_equal = np.abs(K_diff) < 1e-16
        # Kkij_v2k = (Kkij_vinc - Kkij_kees) / Kkij_kees * 100
        # Kkij_k2v = (Kkij_kees - Kkij_vinc) / Kkij_vinc * 100
        # Kkij_log_v2k = np.log(Kkij_v2k)
        # Kkij_log_k2v = np.log(Kkij_k2v)
        
        i = np.linspace(0, N_m, N_m)
        x, y = (i, i) if scale == "lin" else (ac, ac)
        Ks = (K_vinc, K_kees, K_diff)

        plot_1(scale, mg, Ks)
        plot_2(scale, )
      
        # plot_3(cfg, mg, K_vinc)
        # plot_3(cfg, mg, K_kees)

        # from visualization.kernel.kernel import KernelSubplot
        # axis = KernelAxis.Bin if scale == "lin" else KernelAxis.Radius
        # s = KernelSubplot(mg, kernel_vinc.K, axis=axis)
        # p = GridspecPlot([s], add_slider=True)
        # p.render()
