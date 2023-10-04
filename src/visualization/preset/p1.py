import matplotlib.pyplot as plt
import numpy as np

from axis import DiscreteMassAxis, DiscreteTimeAxis, KernelAxisLabelVariant
from config import Config, PATH_TO_DARKMODE
from disk import mass_distribution
from kernel import Kernel
from solver import Solver
from utils.physics import disk_mass_from_distribution  # TODO Remove
from visualization.base import GridspecPlot
from visualization.evolution import EvolutionPlot, MassConservationPlot
from visualization.kernel import KernelSubplot, KernelMassConservationSubplot


def plot_kernel(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    axis: KernelAxisLabelVariant,
    z_limits: tuple[float, float],
):
    if kernel.K.all() >= 0: 
        z_limits = (-z_limits[1], z_limits[1])
    GridspecPlot([
        KernelSubplot(
            cfg, mg, kernel.K,
            title="kernel $K_{kij}$",
            axis_scales=(scale, scale, "lin"),
            z_limits=z_limits,
            symmetrized=True,
            axis_label_variant=axis,
            cmap="bwr",
        )
    ], add_slider=True).render()


def plot_kernel_gain_loss(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    axis: KernelAxisLabelVariant,
    z_limits: tuple[float, float],
):
    cmap = "Reds" if scale == "log" else "bwr"
    GridspecPlot([
        KernelSubplot(
            cfg, mg, kernel.K_gain, 
            title="kernel gain contribution $G_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            symmetrized=True,
            axis_label_variant=axis,
            cmap=cmap,
        ),
        KernelSubplot(
            cfg, mg, -kernel.K_loss,
            title="kernel loss contribution $L_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            ylabel="",
            symmetrized=True,
            axis_label_variant=axis,
            cmap=cmap,
        ),
    ], add_slider=True).render()


def plot_kernel_error(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    axis_label_variant: KernelAxisLabelVariant,
    z_limits: tuple[float, float],
):
    p = GridspecPlot([
        KernelMassConservationSubplot(
            cfg, mg, kernel.K,
            axis_scales=(scale, scale, scale),
            axis_label_variant=axis_label_variant,
            symmetrized=True,
            # z_limits=z_limits, # TODO
        ),
    ])
    p.render()


def integrate(
    cfg: Config,
    kernel: Kernel,
):
    tg = DiscreteTimeAxis(cfg)
    tc = tg.bin_centers
    mc = kernel.mg.bin_centers
    dm = kernel.mg.bin_widths
    K = kernel.K

    n0 = mass_distribution.dirac_delta(cfg)
    solver = Solver(cfg)
    N, f, m2f, dm2f = solver.run(n0, K)
    M = [disk_mass_from_distribution(n, mc, dm) for n in f]

    return tc, f, N, m2f, dm2f, M


def plot_evolution(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    t: np.ndarray,
    N: np.ndarray,
    f: np.ndarray,
    m2f: np.ndarray,
    dm2f: np.ndarray,
):
    # TODO Fix y-limits
    EvolutionPlot(kernel, N, f, m2f, dm2f).render()


def plot_surface(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
    scale: str,
    t: np.ndarray,
    N: np.ndarray,
    f: np.ndarray,
    m2f: np.ndarray,
    dm2f: np.ndarray,
):
    X, Y = mg.bin_centers, t
    X, Y = np.meshgrid(X, Y)
    Z = m2f
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(np.log10(X), np.log10(Y), np.log10(Z))
    ax.set_zlim(-40, -8)
    plt.show()
    plt.close()


def plot_error(
    cfg: Config,
    mg: DiscreteMassAxis,
    kernel: Kernel,
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
    axis_label_variant = KernelAxisLabelVariant.Radius if scale == "log" else KernelAxisLabelVariant.Bin
    z_limits = (1e-20, 1e-7) if scale == "log" else (-1, 1)

    # Plot total kernel     with lin. colorscale.
#    plot_kernel(cfg, mg, kernel, scale, axis_label_variant, z_limits)
    # Plot K_gain & K_loss  with log. colorscale.
    plot_kernel_gain_loss(cfg, mg, kernel, scale, axis_label_variant, z_limits)
    # Plot kernel mass error.
    plot_kernel_error(cfg, mg, kernel, scale, axis_label_variant, z_limits)

#    # Integrate.
#    t, f, N, m2f, dm2f, M = integrate(cfg, kernel)
#
#    # Plot evolution of mass distribution over time.
#    plot_evolution(cfg, mg, kernel, scale, t, f, N, m2f, dm2f)
#    # plot_surface(cfg, mg, kernel, scale, t, f, N, m2f, dm2f)
#
#    # Plot mass error over time.
#    plot_error(cfg, mg, kernel, t, M)
