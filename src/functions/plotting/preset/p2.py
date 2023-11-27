import numpy as np

from config import Config
from models.axis import DiscreteMassAxis
from models.kernel import Kernel
from functions.kernel.mass_conservation import test_mass_conservation
from models.plotting.base import GridspecPlot
from models.plotting.kernel import KernelSubplot, KernelMassConservationSubplot


def plot_sampling_probability_vs_time(
    cfg:            Config, 
    mg:             DiscreteMassAxis, 
    Ps:             np.ndarray,
    symmetrized:    bool                = False
):
    GridspecPlot([
        KernelSubplot(cfg, mg, np.array(Ps), 
            title="Collision Pair Sampling Probability $P_{ij}$",
            symmetrized=symmetrized,
            z_limits=(1e-5, 1),
            cmap="Blues", 
        ),
    ], add_slider=True, slider_label="$i_t$").render()


def plot_sampling_count_vs_time(
    cfg:            Config, 
    mg:             DiscreteMassAxis, 
    Ns:             np.ndarray,
    symmetrized:    bool                = False
):
    GridspecPlot([
        KernelSubplot(cfg, mg, np.array(Ns), 
            title=r"Collision Pair Sampling Count $N_{ij}$",
            axis_scales=("log", "log", "lin"), 
            symmetrized=symmetrized,
            z_limits=(0, np.max(Ns)),  # <- NOTE: Low upper boundary for better visibility.
            cmap="Blues", 
        ),
    ], add_slider=True, slider_label="$i_t$").render()


def plot_kernel_mass_error_vs_time(
    cfg:            Config, 
    mg:             DiscreteMassAxis, 
    Ks:             list[np.ndarray], 
    R_coll:         np.ndarray,
    symmetrized:    bool                = False
):
    Es = [test_mass_conservation(mg, K, R_coll)[0] for K in Ks]
    GridspecPlot([
        # NOTE This should be a `KernelMassConservationSubplot`, but it can't be,
        #      since that class takes `K` as input, not `Es`...
        KernelSubplot(cfg, mg, np.array(Es), 
            title=r"Kernel mass error $\Delta K_{ij}$",
            symmetrized=symmetrized,
            cmap="Reds", 
        ),
    ], add_slider=True, slider_label="$i_t$").render()


def plot_kernel_sampled_vs_unsampled(
    cfg:                Config, 
    mg:                 DiscreteMassAxis, 
    kernel_unsampled:   Kernel, 
    kernel_sampled:     Kernel,
    symmetrized:        bool                = False
):
    K_g, S_g = kernel_unsampled.K_gain, kernel_sampled.K_gain
    K_l, S_l = kernel_unsampled.K_loss, kernel_sampled.K_loss

    cmap, scale, z_limits = "Reds" , "log", (1e-20, 1)
    GridspecPlot([
        KernelSubplot(cfg, mg, K_g,
            title="unsampled kernel gain contr. $G_{kij}$", z_limits=z_limits, 
            axis_scales=(scale, scale, scale), symmetrized=symmetrized, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, -K_l,
            title="unsampled kernel loss contr. $L_{kij}$", z_limits=z_limits, 
            axis_scales=(scale, scale, scale), ylabel="", symmetrized=False, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, S_g,
            title="sampled kernel gain contr. $G_{kij}$", z_limits=z_limits, 
            axis_scales=(scale, scale, scale), symmetrized=False, cmap=cmap,
        ),
        KernelSubplot(cfg, mg, -S_l,
            title="sampled kernel loss contr. $L_{kij}$", z_limits=z_limits, 
            axis_scales=(scale, scale, scale), ylabel="", symmetrized=False, cmap=cmap,
        ),
    ], add_slider=True).render()
