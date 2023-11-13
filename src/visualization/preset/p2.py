import numpy as np

from visualization.base import GridspecPlot
from visualization.kernel import KernelSubplot, KernelMassConservationSubplot
from kernel.mass_conservation import test_mass_conservation


def plot_sampling_probability_vs_time(cfg, mg, Ps):
    GridspecPlot(
        [
            KernelSubplot(cfg, mg, np.array(Ps), 
                cmap="Blues", z_limits=(1e-5, 1),
                title="Collision Pair Sampling Probability $P_{ij}$",
            ),
        ], add_slider=True, slider_label="$i_t$"
    ).render()


def plot_sampling_count_vs_time(cfg, mg, Ns):
    GridspecPlot(
        [
            KernelSubplot(cfg, mg, np.array(Ns), 
                z_limits=(0, np.max(Ns)),  # <- NOTE: Low upper boundary for better visibility.
                axis_scales=("log", "log", "lin"), cmap="Blues", 
                title=r"Collision Pair Sampling Count $N_{ij}$",
            ),
        ], add_slider=True, slider_label="$i_t$"
    ).render()


def plot_kernel_mass_error_vs_time(cfg, mg, Ks, R_coll):
    Es = [test_mass_conservation(mg, K, R_coll)[0] for K in Ks]
    GridspecPlot(
        [
            KernelSubplot(cfg, mg, np.array(Es), 
                cmap="Reds", title=r"Kernel mass error $\Delta K_{ij}$",
            ),
        ], add_slider=True, slider_label="$i_t$"
    ).render()


def plot_kernel_sampled_vs_unsampled(cfg, mg, kernel_unsampled, kernel_sampled):
    K_g, S_g = kernel_unsampled.K_gain, kernel_sampled.K_gain
    K_l, S_l = kernel_unsampled.K_loss, kernel_sampled.K_loss

    cmap, scale, z_limits = "Reds" , "log", (1e-20, 1)
    GridspecPlot([
        KernelSubplot(cfg, mg, K_g,
            title="unsampled kernel gain contr. $G_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits, symmetrized=False, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, -K_l,
            title="unsampled kernel loss contr. $L_{kij}$",
            axis_scales=(scale, scale, scale), ylabel="",
            z_limits=z_limits, symmetrized=False, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, S_g,
            title="sampled kernel gain contr. $G_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits, symmetrized=False, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, -S_l,
            title="sampled kernel loss contr. $L_{kij}$",
            axis_scales=(scale, scale, scale), ylabel="",
            z_limits=z_limits, symmetrized=False, cmap=cmap,
        ),
    ], add_slider=True).render()
