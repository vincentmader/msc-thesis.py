import os
from pathlib import Path
import sys
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_OUTFILES, PATH_TO_FIGURES
from constants import SECONDS_PER_YEAR
from models.solver import SolverV2

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "106")
path_to_figures = Path(PATH_TO_FIGURES, "106")
os.makedirs(path_to_outfiles, exist_ok=True)
os.makedirs(path_to_figures, exist_ok=True)

exports = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k", "S_ij", "P_ij"]
SAMPLING_DENSITIES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]  # NOTE: `1.0` must be first! but then plot is in wrong order.
N_m = 50


def integrate(cfg: Config, title: str):
    path = Path(path_to_outfiles, title)
    if not os.path.exists(path):
        solver = SolverV2(cfg)
        solver.run()
        solver.save(path, exports)

        tc     = solver.tg.bin_centers
        mb_k   = solver.mg.bin_boundaries
        mc_k   = solver.mg.bin_centers
        dm_k   = solver.mg.bin_widths
        n_k    = solver.n_k_vs_t
        N_k    = solver.N_k_vs_t
        M_k    = solver.M_k_vs_t
        dMdt_k = solver.dMdt_vs_t
        P_ij = solver.P_ij_vs_t
        S_ij = solver.S_ij_vs_t
    else:
        tc     = np.loadtxt(Path(path_to_outfiles, title, "tc.txt"))
        mb_k   = np.loadtxt(Path(path_to_outfiles, title, "mb_k.txt"))
        mc_k   = np.loadtxt(Path(path_to_outfiles, title, "mc_k.txt"))
        dm_k   = np.loadtxt(Path(path_to_outfiles, title, "dm_k.txt"))
        n_k    = np.loadtxt(Path(path_to_outfiles, title, "n_k.txt"))
        N_k    = np.loadtxt(Path(path_to_outfiles, title, "N_k.txt"))
        M_k    = np.loadtxt(Path(path_to_outfiles, title, "M_k.txt"))
        dMdt_k = np.loadtxt(Path(path_to_outfiles, title, "dMdt_k.txt"))
        P_ij   = np.loadtxt(Path(path_to_outfiles, title, "P_ij.txt"))
        S_ij   = np.loadtxt(Path(path_to_outfiles, title, "S_ij.txt"))

    return tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k, P_ij, S_ij


if __name__ == "__main__":
    cfgs, titles = [], []
    for rho_sample in SAMPLING_DENSITIES:
        nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
        cfg = Config(
            enable_collision_sampling=True if rho_sample < 1.0 else False,
            nr_of_samples=nr_of_samples,
        )
        cfgs.append(cfg)
        titles.append(f"{rho_sample=}")

    plt.figure(figsize=(10, 11))

    for subplot_idx, rho_sample in tqdm(enumerate(SAMPLING_DENSITIES)):
        plt.subplot(4, 2, subplot_idx+1)

        cfg = cfgs[subplot_idx]
        title = titles[subplot_idx]
        tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k, P_ij, S_ij = integrate(cfg, title)
        M_k_complete = np.loadtxt(Path(path_to_outfiles, "rho_sample=1.0", "M_k.txt"))

        X_MIN = tc[0] / SECONDS_PER_YEAR
        X_MAX = tc[-1] / SECONDS_PER_YEAR
        
        err = (M_k.sum(axis=1) - M_k[0].sum(axis=0)) / M_k[0].sum(axis=0)
        plt.loglog(tc / SECONDS_PER_YEAR, err, 'g', 
            # label=r"$\Delta_{stab}(t)=\sum_i \frac{\rho^d_i(t) - \rho^d_i(t=0)}{\rho^d_i(t=0)}$")
            label=r"$\Delta_{stab}(t)$")
        plt.loglog(tc / SECONDS_PER_YEAR, -err, 'g--')

        if rho_sample != 1.0:
            # err = (((M_k - M_k_complete) / M_k_complete)**2).sum(axis=1)**.5
            err = (((M_k - M_k_complete) / M_k_complete)**2).sum(axis=1)**.5
            # err = np.array([np.correlate(M_k, M_k_complete) for M_k, M_k_complete in zip(M_k, M_k_complete)])

            # a = err[err.nan()]
            for i in err:
                print(i)

            plt.loglog(tc / SECONDS_PER_YEAR, err, 'b', 
                label=r"$\Delta_{acc}(t)$")
            plt.loglog(tc / SECONDS_PER_YEAR, -err, 'b--')

        plt.title(r"$\rho_{sample} = "+ f"{rho_sample}" + "$")
        if subplot_idx + 2 >= len(SAMPLING_DENSITIES):
            plt.xlabel("Time $t$ [y]")
        if subplot_idx % 2 == 0:
            # plt.ylabel(r"error $\Delta_{acc}$ and $\Delta_{stab}$")
            plt.ylabel("Error (Dimensionless)")

        # plt.plot([X_MIN, X_MAX], [1, 1], color="black", label="1")

        if subplot_idx == 0:
            plt.legend(loc="upper left")
        plt.xlim(X_MIN, X_MAX)
        plt.ylim(1e-17, 1e14)
        plt.grid(True)
        plt.tight_layout()

    plt.savefig(Path(path_to_figures, "accuracy_and_stability_vs_sampling_density.pdf"))
    plt.show()
    plt.close()
