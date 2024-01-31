import os, sys
from pathlib import Path
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
times_of_interest = list(range(46, 190, 1))
# times_of_interest = [0, 49, 99, 149, 159, 169, 179, 189, 199]
N_m = 50


def plot_mass_distribution_vs_m_vs_t(title, mb_k, mc_k, M_k, should_show_plot=False):

    plt.figure(figsize=(11, 7))
    colors = plt.cm.Blues(np.linspace(-0.1, 1.1, len(times_of_interest)))
    x = mc_k

    for ii, i_t in enumerate(times_of_interest):
        y = M_k[i_t]
        plt.loglog(
            x, y,
            color=colors[ii]
        )

    plt.title(title)
    plt.xlim(mb_k[0], mb_k[-1])
    plt.ylim(1e-20, 1e-5)
    plt.grid(True)

    name = f"mass_distr {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()

def plot_mass_error_vs_t(title, tc, M_k, should_show_plot=False):
    plt.figure(figsize=(11, 4))

    M_t = np.sum(M_k, axis=1)
    M_0 = M_t[0]

    x = tc / SECONDS_PER_YEAR
    y = np.array([(M_t - M_0) / M_0 for M_t in M_t])
    plt.loglog(x, y, color="red", label=r"$\Delta\rho^{rel}_d > 0$")
    plt.loglog(x, -y, color="blue", label=r"$\Delta\rho^{rel}_d < 0$")

    plt.ylim(1e-16, 1e-12)
    plt.xlabel("time [y]")
    plt.ylabel(r"relative density error $\frac{ \rho_d(t) - \rho_d(t=0) }{ \rho_d(t=0) }$")
    plt.grid(True)

    plt.legend(loc="upper left")

    name = f"mass_error {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()


def plot_rms_mass_deriv_vs_t(title, tc, dMdt_k, should_show_plot=False):
    plt.figure(figsize=(11, 4))

    # M_t = np.sum(M_k, axis=1)
    # M_0 = M_t[0]

    y = (dMdt_k**2).sum(axis=1)**.5

    x = tc / SECONDS_PER_YEAR
    # y = np.array([(M_t - M_0) / M_0 for M_t in M_t])
    plt.loglog(x, y, color="red", label=r"$\sqrt{ \sum_i \left( \frac{\Delta n_i}{\Delta t}\cdot m_i\cdot \Delta m_i \right)^2 }$")

    plt.ylim(1e-32, 1e-12)
    plt.xlabel("time [y]")
    # plt.ylabel(r"relative density error $\frac{ \rho_d(t) - \rho_d(t=0) }{ \rho_d(t=0) }$")
    plt.grid(True)

    plt.legend(loc="upper left")

    name = f"mass_deriv {title}.pdf"
    path = Path(path_to_figures, name)
    plt.savefig(path)
    if should_show_plot:
        plt.show()
    plt.close()
    

def foo(cfg, title, should_show_plot=False):
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

    plot_mass_distribution_vs_m_vs_t(title, mb_k, mc_k, M_k, should_show_plot=should_show_plot)
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
        foo(cfg, title, should_show_plot=should_show_plot)

if __name__ == "__main__":
    main(
        should_show_plot=True,
    )
