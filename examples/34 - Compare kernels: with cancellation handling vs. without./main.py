import os, sys
import numpy as np
from pathlib import Path
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import AxisLabelVariant
    from config import Config, PATH_TO_FIGURES
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel import KernelMassConservationSubplot, KernelSubplot
except ModuleNotFoundError as e:
    raise e

kwargs = {
    "mass_max_value": 1e10, 
    "mass_resolution": 100
}
cfg_1 = Config(enable_cancellation_handling=False, **kwargs)
cfg_2 = Config(enable_cancellation_handling=True, **kwargs)
kernel_1, kernel_2 = Kernel(cfg_1), Kernel(cfg_2)
K_1, K_2 = kernel_1.K, kernel_2.K
K_diff = K_1 - K_2
K_equal = K_diff < 1e-14

mg = kernel_1.mg
mc = mg.bin_centers
ac = mg.particle_radii

R_coll = np.ones(shape=[mg.N]*2)


def plot_1():

    # Plot kernels `K_canc` `K_nocanc` side by side, & plot difference.
    GridspecPlot([
        KernelSubplot(
            cfg_1, mg, K_1,
            title="$K_{kij}^{canc}$",
            axis_scales=("log", "log", "lin"),
            symmetrized=True,
            cmap="bwr",
        ),
        KernelSubplot(
            cfg_2, mg, K_2,
            title="$K_{kij}^{nocanc}$ with canc. handling",
            axis_scales=("log", "log", "lin"),
            symmetrized=True,
            cmap="bwr",
        ),
        KernelSubplot(
            cfg_1, mg, np.abs(K_diff), # TODO This is the wrong `cfg`.
            title="abs($K_{kij}^{canc}-K_{kij}^{nocanc}$)",
            axis_scales=("log", "log", "log"),
            symmetrized=True,
    )], add_slider=True).render()


def plot_2():

    # Plot kernel errors of `K_canc` & `K_nocanc` side by side.
    s1 = KernelMassConservationSubplot(
        cfg_1, mg, K_1, R_coll,
        axis_label_variant=AxisLabelVariant.Radius,
        title=r"$\Delta K^{canc}_{ij}$",
    )
    s2 = KernelMassConservationSubplot(
        cfg_2, mg, K_2, R_coll,
        axis_label_variant=AxisLabelVariant.Radius,
        title=r"$\Delta K^{nocanc}_{ij}$",
    )
    GridspecPlot([s1, s2]).render()

    # Plot kernel errors of `K_canc` & `K_nocanc` separately, one after the other.
    for label, subplot in zip(["canc", "nocanc"], [s1, s2]):
        path = Path(PATH_TO_FIGURES, "34")
        os.makedirs(path, exist_ok=True)
        path = Path(path, f"{label}.pdf")
        GridspecPlot([subplot]).render(
            show_plot=False, save_plot=True, path_to_outfile=path
        )


def plot_3():
    for (enable_coagulation, enable_fragmentation) in [(True, False), (False, True)]:

        kernel_superscript = ""
        kernel_superscript = "coag" if (enable_coagulation and not enable_fragmentation) else kernel_superscript
        kernel_superscript = "frag" if (not enable_coagulation and enable_fragmentation) else kernel_superscript

        for enable_cancellation_handling in [True, False]:

            kernel_subscript = "nocanc" if enable_cancellation_handling else "canc"

            cfg = Config(
                enable_coagulation=enable_coagulation,
                enable_fragmentation=enable_fragmentation,
                enable_cancellation_handling=enable_cancellation_handling,
            )
            kernel = Kernel(cfg)
            mg, K = kernel.mg, kernel.K
            # E_ij, E_tot = test_mass_conservation(mg, K)

            title = r"K^{" + kernel_superscript + "}_{" + kernel_subscript + "}"
            filename = title.replace("{", "").replace("}", "").replace("^", "_")
            filename = f"{filename}.pdf"
            path_to_outfile = Path(PATH_TO_FIGURES, "34", filename)
            title = f"${title}$"

            GridspecPlot([
                KernelMassConservationSubplot(
                    cfg, mg, K, R_coll,
                    title=title,
                ),
            ]) .render(save_plot=True, path_to_outfile=path_to_outfile)


if __name__ == "__main__":
    plot_1()
    plot_2()
    plot_3()
