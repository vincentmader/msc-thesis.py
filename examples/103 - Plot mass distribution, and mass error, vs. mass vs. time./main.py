import os, sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from models.kernel import Kernel
from models.solver import SolverV2
from constants import SECONDS_PER_YEAR

path_to_figures  = Path(PATH_TO_FIGURES, "103")
path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "103")
os.makedirs(path_to_figures, exist_ok=True)

exports = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k"]
N_m = 50


def plot_mass_distribution_vs_m_vs_t(title, tc, mb_k, mc_k, M_k, should_show_plot=False, enable_coagulation=None, enable_fragmentation=None):
    if enable_coagulation is False and enable_fragmentation is True:
        times_of_interest = [0] + list(range(50, 190, 5))
        # times_of_interest = list(range(150, 190, 2))
        # times_of_interest = [0, 49, 99, 149, 159, 169, 179, 189, 199]
    else:
        times_of_interest = [0] + list(range(100, 190, 5))

    fig = plt.figure(figsize=(10, 5))
    colors = plt.cm.YlGnBu(np.linspace(0, 1, len(times_of_interest)))
    x = mc_k

    for ii, i_t in enumerate(times_of_interest):
        y = M_k[i_t]
        plt.loglog(
            x, y,
            color=colors[ii]
        )

    # plt.title(title)
    plt.xlim(mb_k[0], mb_k[-1])
    plt.ylim(1e-16, 1e-8)
    plt.xlabel("Dust Particle Mass $m_i$ [kg]")
    plt.ylabel(r"Dust Density $\rho^d_i$ [kg m$^{-3}$]")
    plt.grid(True)

    plt.subplots_adjust(top=0.85)
    # cax = plt.axes([0.91, 0.1, 0.021, 0.78])
    cax = fig.add_axes([0.125, 0.9, 0.775, 0.05])
    cb = mpl.colorbar.ColorbarBase(cax, 
        orientation="horizontal", 
        cmap="YlGnBu",
        norm=mpl.colors.LogNorm(
            tc[times_of_interest[0]] / SECONDS_PER_YEAR, 
            tc[times_of_interest[-1]] / SECONDS_PER_YEAR,
        ),
    )
    cb.ax.set_title("Time $t$ [y]")

    name = f"mass_distr {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()

def plot_mass_error_vs_t(title, tc, M_k, should_show_plot=False):
    plt.figure(figsize=(10, 4))

    M_t = np.sum(M_k, axis=1)
    M_0 = M_t[0]

    x = tc / SECONDS_PER_YEAR
    y = np.array([(M_t - M_0) / M_0 for M_t in M_t])
    plt.loglog(x, y, color="red", label=r"$\Delta_{stab} > 0$")
    plt.loglog(x, -y, color="blue", label=r"$\Delta_{stab} < 0$")

    plt.ylim(1e-17, 1e-7)
    plt.xlabel("Time [y]")
    plt.ylabel(r"Dimensionless Mass Error $\Delta_{stab}(t)$")
    plt.grid(True)

    plt.legend(loc="upper left")

    name = f"mass_error {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()


def plot_rms_mass_deriv_vs_t(title, tc, dMdt_k, should_show_plot=False):
    plt.figure(figsize=(10, 4))

    # M_t = np.sum(M_k, axis=1)
    # M_0 = M_t[0]

    y = (dMdt_k**2).sum(axis=1)**.5

    x = tc / SECONDS_PER_YEAR
    # y = np.array([(M_t - M_0) / M_0 for M_t in M_t])
    plt.loglog(
        x, y, color="red", 
        label=r"$\left(\frac{\partial \rho}{\partial t}\right)_{RMS} = \sqrt{ \sum_i \left( \frac{\Delta n_i}{\Delta t}\cdot m_i\cdot \Delta m_i \right)^2 }$"
    )

    plt.ylim(1e-32, 1e-12)
    plt.xlabel("Time [y]")
    plt.ylabel("RMS of Temporal Density Derivative [kg m$^{-3}$ s$^{-1}$]")
    # plt.ylabel(r"relative density error $\frac{ \rho_d(t) - \rho_d(t=0) }{ \rho_d(t=0) }$")
    plt.grid(True)

    plt.legend(loc="upper left")

    name = f"mass_deriv {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()
    

def foo(cfg, title, should_show_plot=False, enable_coagulation=None, enable_fragmentation=None):
    path = Path(path_to_outfiles, title)
    if not os.path.exists(path):
        solver = SolverV2(cfg)
        solver.run()
        solver.save(path, exports=exports)
        tc     = solver.tg.bin_centers
        mb_k   = solver.mg.bin_boundaries
        mc_k   = solver.mg.bin_centers
        dm_k   = solver.mg.bin_widths
        n_k    = solver.n_k_vs_t
        N_k    = solver.N_k_vs_t
        M_k    = solver.M_k_vs_t
        dMdt_k = solver.dMdt_vs_t
    else:
        tc     = np.loadtxt(Path(path_to_outfiles, title, "tc.txt"))
        mb_k   = np.loadtxt(Path(path_to_outfiles, title, "mb_k.txt"))
        mc_k   = np.loadtxt(Path(path_to_outfiles, title, "mc_k.txt"))
        dm_k   = np.loadtxt(Path(path_to_outfiles, title, "dm_k.txt"))
        n_k    = np.loadtxt(Path(path_to_outfiles, title, "n_k.txt"))
        N_k    = np.loadtxt(Path(path_to_outfiles, title, "N_k.txt"))
        M_k    = np.loadtxt(Path(path_to_outfiles, title, "M_k.txt"))
        dMdt_k = np.loadtxt(Path(path_to_outfiles, title, "dMdt_k.txt"))

    plot_mass_distribution_vs_m_vs_t(title, tc, mb_k, mc_k, M_k, should_show_plot=should_show_plot,
                                     enable_coagulation=enable_coagulation, enable_fragmentation=enable_fragmentation)
    plot_mass_error_vs_t(title, tc, M_k, should_show_plot=should_show_plot)
    plot_rms_mass_deriv_vs_t(title, tc, dMdt_k, should_show_plot=should_show_plot)

def main(
    should_show_plot=False,
):
    cfgs = []
    for (enable_coagulation, enable_fragmentation) in [(True, True), (True, False), (False, True)]:
        for m_0 in [0, int(N_m / 2)]:
            if enable_coagulation == True and m_0 != 0: 
                continue
            cfg = Config(
                enable_collision_sampling=False,
                time_max_value=1e11,
                enable_coagulation=enable_coagulation,
                enable_fragmentation=enable_fragmentation,
                mass_resolution=N_m,
                initial_mass_bin=m_0,
            )
            cfgs.append(cfg)

    for cfg in cfgs:
        m_0 = cfg.initial_mass_bin
        title = f"coag={cfg.enable_coagulation} frag={cfg.enable_fragmentation} m0={m_0}"
        foo(cfg, title, should_show_plot=should_show_plot, 
            enable_coagulation=cfg.enable_coagulation, enable_fragmentation=cfg.enable_fragmentation)

if __name__ == "__main__":
    main(
        should_show_plot=False,
    )
