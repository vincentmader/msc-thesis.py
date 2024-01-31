import os, sys
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from constants import SECONDS_PER_YEAR
from models.axis import DiscreteMassAxis
from models.solver import SolverV2
from functions.utils.dates import format_seconds_as_years


path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "109")
path_to_figures  = Path(PATH_TO_FIGURES, "109")
os.makedirs(path_to_figures, exist_ok=True)
exports = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k", "S_ij", "P_ij", "X_ij"]
# times_of_interest = [int(i_t) for i_t in np.linspace(0, 199, 12)]
# times_of_interest = [1, 109, 129, 149, 159, 169, 174, 179, 184, 189, 194, 199]


def foo(cfg: Config, title: str, label: str, color=None):

    mg = DiscreteMassAxis(cfg)

    path = Path(path_to_outfiles, title)
    if not os.path.exists(path):
        solver = SolverV2(cfg)
        solver.run()
        solver.save(path, exports)

        tc     = solver.tg.bin_centers
        X_ij   = solver.X_ij_vs_t
    else:
        tc     = np.loadtxt(Path(path_to_outfiles, title, "tc.txt"))
        X_ij   = np.loadtxt(Path(path_to_outfiles, title, "X_ij.txt"))

    x = tc / SECONDS_PER_YEAR
    y = X_ij
    # plt.semilogx(x, y, label=label, color=color)
    plt.loglog(x, y, label=label, color=color)


def main(should_show=False):
#     cfgs = []

    SAMPLING_DENSITIES = [0.2, 0.4, 0.6, 0.8, 1.0]
    # colors = ["cyan", "green", "orange", "red", "black"]
    TITLES = [
        "full model",
        "pure coagulation",
        "only fragmentation",
    ]
    N_m = 100

    for i, (enable_coagulation, enable_fragmentation) in enumerate([(True, True), (True, False), (False, True)]):

        plt.figure(figsize=(10, 4))

        for ii, rho_sample in enumerate(SAMPLING_DENSITIES):
            color = "k" if rho_sample == 1 else None

            nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
            cfg = Config(
                mass_resolution=N_m,
                nr_of_samples=nr_of_samples,
                enable_collision_sampling=(rho_sample != 1),
                enable_coagulation=enable_coagulation,
                enable_fragmentation=enable_fragmentation,
            )
            title = f"coag={enable_coagulation} frag={enable_fragmentation} {rho_sample=}"
            label = r"$\rho_{sample} = $" + f"{rho_sample}"
            foo(cfg, title, color=color, label=label)

        plt.grid()
        plt.xlabel("time $t$ [y]")
        plt.ylabel("share of non-zero kernel entries")
        plt.ylim(1e-6, 1)
#         plt.ylabel(r"temporal derivative $\sqrt{\sum_i \left(\frac{\Delta \rho_i}{\Delta t}\right)^2}$")
        plt.legend(loc="lower left")
        plt.title(TITLES[i])

        plt.savefig(Path(path_to_figures, f"percentage non-zero kernel, coag={enable_coagulation}, frag={enable_fragmentation}.pdf"))
        if should_show:
            plt.show()
        plt.close()


if __name__ == "__main__":
    main(should_show=True)
