import os, sys
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "lib"))
    from coag_py.coag.coag import create_coag_kernel
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, AxisLabelVariant
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelSubplot, KernelMassConservationSubplot
except ModuleNotFoundError as e:
    raise e


def plot_1(
    cfg: Config,
    scale: str,
    mg: DiscreteMassAxis,
    Ks: tuple[np.ndarray, np.ndarray, np.ndarray],
):
    K_vinc, K_kees, K_compare = Ks
    axis_label_variant = AxisLabelVariant.Bin if scale == "lin" else AxisLabelVariant.Radius

    GridspecPlot([
        KernelSubplot(
            cfg, mg, K_vinc,
            title="$K_{kij}^{vinc}$",
            axis_label_variant=axis_label_variant,
            symmetrized=True,
            axis_scales=(scale, scale, "lin"),
            z_limits=(-1, 1),
            cmap="bwr",
        ),
        KernelSubplot(
            cfg, mg, K_kees,
            title="$K_{kij}^{kees}$",
            axis_label_variant=axis_label_variant,
            symmetrized=True,
            axis_scales=(scale, scale, "lin"),
            z_limits=(-1, 1),
            cmap="bwr",
        ),
        KernelSubplot(
            cfg, mg, K_compare,
            title=r"$\Delta K_{kij}=K_{kij}^{kees}-K_{kij}^{vinc}$",
            axis_label_variant=axis_label_variant,
            symmetrized=True,
            axis_scales=(scale, scale, "lin"),
            z_limits=(-1, 1),
            cmap="bwr",
        ),
    ], add_slider=True).render()


def plot_2(
    cfg: Config,
    scale: str,
    mg: DiscreteMassAxis,
    Ks: tuple[np.ndarray, np.ndarray],
):
    K_vinc, K_kees = Ks
    axis_label_variant = AxisLabelVariant.Bin if scale == "lin" else AxisLabelVariant.Radius

    s1 = KernelMassConservationSubplot(
        cfg, mg, K_vinc, axis_label_variant=axis_label_variant, axis_scales=(scale, scale, "log"), 
    )
    s2 = KernelMassConservationSubplot(
        cfg, mg, K_kees, axis_label_variant=axis_label_variant, axis_scales=(scale, scale, "log"),
    )
    p = GridspecPlot([s1, s2])
    p.render()


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
            enable_cancellation_handling=True,
        )
    
        kernel_vinc = Kernel(cfg)
        K_vinc, mg = kernel_vinc.K, kernel_vinc.mg
        mc, ac = mg.bin_centers, mg.particle_radii
    
        R_coll = np.ones(shape=[N_m] * 2) # TODO Redefine `R` (and `K` ?)
        K_kees = create_coag_kernel(mc, R_coll)
        
        K_vinc_sym = np.array([0.5 * (z + z.T) for z in K_vinc])
        K_diff = K_kees - K_vinc_sym
        
        i = np.linspace(0, N_m, N_m)
        x, y = (i, i) if scale == "lin" else (ac, ac)

        plot_1(cfg, scale, mg, (K_vinc, K_kees, K_diff))  # TODO Fix this!!!
        plot_2(cfg, scale, mg, (K_vinc, K_kees))
