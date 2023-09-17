import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_radius_from_mass, particle_mass_from_radius
    from kees_kernel import create_coag_kernel
    from kernel import Kernel
    from kernel.mass_conservation import test_mass_conservation
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
    from visualization.kernel.mass_conservation import KernelMassConservationSubplot
    from visualization.v1.mass_conservation import KernelMassConservationPlot
except ModuleNotFoundError as e:
    raise e


def plot_1():
    pass

def plot_2():
    pass

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
        K_vinc = kernel_vinc.K
    
        mg = kernel_vinc.mg
        mc = mg.grid_cell_centers
        rho_s = cfg.dust_particle_density
        ac = particle_radius_from_mass(mc, rho_s)
    
        R_coll = np.ones(shape=[N_m] * 2) # TODO Redefine `R` (and `K` ?)
        K_kees = create_coag_kernel(mc, R_coll)
        
        N_m = mg.N
        i = np.linspace(0, N_m, N_m)
        
        K_vinc_sym = np.array([0.5 * (z + z.T) for z in K_vinc])
        K_diff = K_vinc_sym - K_kees
        K_equal = np.abs(K_diff) < 1e-16
        # Kkij_v2k = (Kkij_vinc - Kkij_kees) / Kkij_kees * 100
        # Kkij_k2v = (Kkij_kees - Kkij_vinc) / Kkij_vinc * 100
        # Kkij_log_v2k = np.log(Kkij_v2k)
        # Kkij_log_k2v = np.log(Kkij_k2v)
        
        s1 = PcolorMatrixSubplot(
            i, i, K_vinc,
            title="$K_{kij}^{vinc}$",
            xlabel="bin index $j$",
            ylabel="bin index $i$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            z_limits=(-1, 1),
        )
        s2 = PcolorMatrixSubplot(
            i, i, K_kees,
            title="$K_{kij}^{kees}$",
            xlabel="bin index $j$",
            scales=("lin", "lin", "lin"),
            z_limits=(-1, 1),
        )
        s3 = PcolorMatrixSubplot(
            i, i, K_equal,
            title="$K_{kij}^{vinc}=K_{kij}^{kees}$",
            xlabel="bin index $j$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
        )
        p = GridspecPlot([s1, s2, s3], add_slider=True)
        p.render()
      
        sum_ij = test_mass_conservation(cfg, mg, kernel_vinc.K)
        def custom_format_coord(x, y):
            x, y = particle_mass_from_radius(x, rho_s), particle_mass_from_radius(y, rho_s)
            x, y = mg.index_from_value(x), mg.index_from_value(y)
            i, j = int(x), int(y)  # TODO Index Convention? (irrelevant due to symmetry)
            text = f"sum_k K_kij = {sum_ij[i, j]:.2}, {i = }, {j = }"
            return text
    
        s1 = KernelMassConservationSubplot(
            kernel_vinc,
            title=r"kernel mass error $\sum_k m_k K_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            scales=(scale, scale, "log"),
        )
        z = test_mass_conservation(cfg, mg, K_kees)
        x, y = (i, i) if scale == "lin" else (ac, ac)
        s2 = PcolorMatrixSubplot(
            x, y, z,
            title=r"kernel mass error $\sum_k m_k K_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            scales=(scale, scale, "log"),
        )
        p = GridspecPlot([s1, s2])
        p.axes[1].format_coord = custom_format_coord
        p.render()

        plot_3(cfg, mg, K_vinc)
        plot_3(cfg, mg, K_kees)
