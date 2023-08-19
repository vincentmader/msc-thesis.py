import os
import sys
import matplotlib.pyplot as plt
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config, PATH_TO_DARKMODE, PATH_TO_FIGURES
    from axis import DiscreteMassAxis
except ModuleNotFoundError as e:
    raise e

FIGSIZE = (5, 5)


def plot_1(i, m, mg, scale):
    if scale == "lin":
        plt.plot(i, m)
    else:
        plt.semilogy(i, m)
    plt.title("dust particle mass vs. bin index")
    plt.xlabel("bin index $i$")
    plt.ylabel("dust particle mass $m_i$")
    plt.xlim(i[0], i[-1])
    plt.ylim(m[0], m[-1])
    plt.grid(True)


def plot_2(i, m, mg, scale):
    if scale == "lin":
        plt.plot(m, i)
    else:
        plt.semilogx(m, i)
    plt.title("bin index vs. dust particle mass")
    plt.xlabel("dust particle mass $m_i$")
    plt.ylabel("bin index $i$")
    plt.xlim(m[0], m[-1])
    plt.ylim(i[0], i[-1])
    plt.grid(True)


def plot_3(i, _, mg, _scale):
    y = [mg.index_from_value(mg.value_from_index(i)) for i in i]
    plt.plot(i, y)
    plt.title(r"$i\to m\to i$ conversion")
    plt.xlabel("bin index $i$")
    plt.ylabel("bin index $i$")
    plt.axis('scaled')
    plt.xlim(i[0], i[-1])
    plt.ylim(y[0], y[-1])
    plt.grid(True)


def plot_4(_, m, mg, scale):
    y = [mg.value_from_index(mg.index_from_value(m)) for m in m]
    if scale == "lin":
        plt.plot(m, y)
    else:
        plt.loglog(m, y)
    plt.title(r"$m\to i\to m$ conversion")
    plt.xlabel("dust particle mass $m_i$")
    plt.ylabel("dust particle mass $m_i$")
    plt.axis('scaled')
    plt.xlim(m[0], m[-1])
    plt.ylim(y[0], y[-1])
    plt.grid(True)


def plot(i, m, mg, scale, plot_separately=False, show_plot=False):
    os.makedirs("../../figures/11", exist_ok=True)
    plots = [plot_1, plot_2, plot_3, plot_4]
    if plot_separately:
        for idx, plot in enumerate(plots):
            plt.figure(figsize=FIGSIZE)
            plot(i, m, mg, scale)
            plt.tight_layout()
            filename = f"discrete-mass-axis_{scale}-scale-{idx}.pdf"
            path = os.path.join(PATH_TO_FIGURES, "11", filename)
            plt.savefig(path)
            plt.gcf().subplots_adjust(left=0.15)
            if show_plot:
                plt.show()
            plt.close()
    else:
        plt.figure(figsize=(9, 9))
        for idx, plot in enumerate(plots):
            plt.subplot(2, 2, idx + 1)
            plot(i, m, mg, scale)
        plt.tight_layout()
        filename = f"discrete-mass-axis_{scale}-scale.pdf"
        path = os.path.join(PATH_TO_FIGURES, "11", filename)
        plt.savefig(path)
        if show_plot:
            plt.show()
        plt.close()


def main():
    cfg = Config()
    if cfg.mpl_dark_mode:
        plt.style.use(PATH_TO_DARKMODE)
    for scale in ["lin", "log"]:
        cfg.mass_axis_scale = scale

        mg = DiscreteMassAxis(cfg)
        i = np.arange(0, mg.N, 1)
        m = mg.grid_cell_centers
        # ^ NOTE: Used bounds (not centers) here,
        #   since centers are indexed by x.5 values.
        plot(i, m, mg, scale, plot_separately=False, show_plot=True)
        plot(i, m, mg, scale, plot_separately=True, show_plot=False)
