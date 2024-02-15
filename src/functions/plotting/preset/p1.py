from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from config import Config, PATH_TO_DARKMODE
from models.axis import DiscreteMassAxis, DiscreteTimeAxis, AxisLabelVariant
from models.kernel import Kernel
from models.plotting.base import GridspecPlot
from models.plotting.evolution import EvolutionPlot, MassConservationPlot
from models.plotting.kernel import KernelSubplot, KernelMassConservationSubplot
from functions.disk import mass_distribution
from functions.utils.physics import disk_mass_from_distribution  # TODO Remove


def plot_kernel(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    z_limits: tuple[float, float],
    axis_label_variant: Optional[AxisLabelVariant] = AxisLabelVariant.Radius,
    symmetrized=False,
):
    if not kernel.K.all() >= 0: 
        z_limits = (-z_limits[1], z_limits[1])
    GridspecPlot([
        KernelSubplot(
            cfg, mg, kernel.K,
            title="kernel $K_{kij}$",
            axis_scales=(scale, scale, "lin"),
            z_limits=z_limits,
            symmetrized=symmetrized,
            axis_label_variant=axis_label_variant,
            cmap="bwr",
        )
    ], add_slider=True).render()


def plot_kernel_gain_loss(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    z_limits: tuple[float, float],
    axis_label_variant: Optional[AxisLabelVariant] = AxisLabelVariant.Radius,
    symmetrized=False,
):
    cmap = "Reds" if scale == "log" else "bwr"
    GridspecPlot([
        KernelSubplot(
            cfg, mg, kernel.K_gain, 
            title=r"kernel gain $G_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            symmetrized=symmetrized,
            axis_label_variant=axis_label_variant,
            cmap=cmap,
        ),
        KernelSubplot(
            cfg, mg, -kernel.K_loss,
            title=r"kernel loss $L_{kij}$, $K_{kij}\cdot R_{ij}$ [$kg m^3 s^{-1}$]",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            ylabel="",
            symmetrized=symmetrized,
            axis_label_variant=axis_label_variant,
            cmap=cmap,
        ),
    ], add_slider=True).render()


def plot_kernel_error(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    R: np.ndarray,
    scale: str,
    z_limits: tuple[float, float],
    axis_label_variant: Optional[AxisLabelVariant] = AxisLabelVariant.Radius,
    symmetrized=False,
):
    GridspecPlot([
        KernelMassConservationSubplot(
            cfg, mg, kernel.K, R,
            axis_scales=(scale, scale, scale),
            axis_label_variant=axis_label_variant,
            symmetrized=symmetrized,
            # z_limits=z_limits, # TODO
        ),
    ]).render()


def integrate(
    cfg: Config,
    kernel: Kernel,
):
    tg = DiscreteTimeAxis(cfg)
    tc = tg.bin_centers
    mc = kernel.mg.bin_centers
    dm = kernel.mg.bin_widths

    from models.solver import SolverV2
    solver = SolverV2(cfg)
    solver.run()

    n_k_vs_t  = solver.n_k_vs_t
    N_k_vs_t  = solver.N_k_vs_t
    M_k_vs_t  = solver.M_k_vs_t
    dMdt_vs_t = solver.dMdt_vs_t

    M = np.array([disk_mass_from_distribution(n, mc, dm) for n in n_k_vs_t])

    return tc, n_k_vs_t, N_k_vs_t, M_k_vs_t, dMdt_vs_t, M


def plot_evolution(
    cfg: Config,
    mg: DiscreteMassAxis,
    scale: str,
    t: np.ndarray,
    f: np.ndarray,
    N: np.ndarray,
    m2f: np.ndarray,
    dm2fdt: np.ndarray,
):
    # TODO Fix y-limits?
    EvolutionPlot(cfg, mg, N=N, f=f, m2f=m2f, dm2fdt=dm2fdt).render()


def plot_surface(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    t: np.ndarray,
    N: np.ndarray,
    f: np.ndarray,
    m2f: np.ndarray,
    dm2fdt: np.ndarray,
):
    X, Y = mg.bin_centers, t
    X, Y = np.meshgrid(X, Y)
    Z = m2f
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(np.log10(X), np.log10(Y), np.log10(Z))
    ax.set_zlim(-20, -8)
    plt.show()
    plt.close()


def plot_error(
    cfg: Config,
    mg: DiscreteMassAxis,
    t: np.ndarray,
    M: np.ndarray,
):
    # TODO Fix y-limits
    MassConservationPlot(cfg, t, M).render()


def main(cfg):
    if cfg.mpl_dark_mode:
        plt.style.use(PATH_TO_DARKMODE)
    mg, kernel = DiscreteMassAxis(cfg), Kernel(cfg)
    scale = mg.scale
    axis_label_variant = AxisLabelVariant.Radius if scale == "log" else AxisLabelVariant.Bin
    z_limits = (1e-20, 1e+10) if scale == "log" else (-1, 1)

    # Plot total kernel     with lin. colorscale.
#    plot_kernel(cfg, mg, kernel, scale, axis_label_variant, z_limits)

    # Plot K_gain & K_loss  with log. colorscale.
    plot_kernel_gain_loss(cfg, mg, kernel, scale, z_limits, axis_label_variant=axis_label_variant)

    # Plot kernel mass error.
    R = kernel.R_coag + kernel.R_frag
    plot_kernel_error(cfg, mg, kernel, R, scale, z_limits, axis_label_variant=axis_label_variant)

    # Integrate.
    t, f, N, m2f, dm2fdt, M = integrate(cfg, kernel)

    # Plot evolution of mass distribution over time.
    plot_evolution(cfg, mg, scale, t, f, N, m2f, dm2fdt)
    # plot_surface(cfg, mg, kernel, scale, t, f, N, m2f, dm2fdt)

    # Plot mass error over time.
    plot_error(cfg, mg, t, M)
