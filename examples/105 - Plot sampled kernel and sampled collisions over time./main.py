import os, sys
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from models.axis import DiscreteMassAxis
from models.solver import SolverV2
from functions.utils.dates import format_seconds_as_years


path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "105")
path_to_figures  = Path(PATH_TO_FIGURES, "105")
os.makedirs(path_to_figures, exist_ok=True)
exports = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k", "S_ij", "P_ij"]
times_of_interest = [int(i_t) for i_t in np.linspace(0, 199, 12)]
times_of_interest = [1, 109, 129, 149, 159, 169, 174, 179, 184, 189, 194, 199]


def main(cfg: Config, title: str, should_show_plot=False):

    mg = DiscreteMassAxis(cfg)
    ac = mg.particle_radii

    N_m = cfg.mass_resolution
    N_t = cfg.time_resolution

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

        P_ij = P_ij.reshape([N_t, N_m, N_m])
        S_ij = S_ij.reshape([N_t, N_m, N_m])

    for symbol, Z in {"S_ij": S_ij, "P_ij": P_ij}.items():
        if symbol == "P_ij":
            cmap = "Blues"
            vmin, vmax = 1e-10, 1 / len(N_k**2)
            norm = LogNorm(vmin=vmin, vmax=vmax)
        elif symbol == "S_ij":
            cmap = "Greens"
            vmin, vmax = 0, Z.max()
            norm = Normalize(vmin=vmin, vmax=vmax)
        else:
            raise Exception()

        fig = plt.figure(figsize=(10, 11))

        for subplot_idx, i_t in enumerate(times_of_interest):
            Z_t = Z[i_t]
            for i in range(N_m):
                for j in range(N_m):
                    if i < j:
                        Z_t[i, j] = Z_t[j, i]
            # Z_t = 0.5 * (Z_t + Z_t.T)

            ax = plt.subplot(4, 3, subplot_idx + 1)

            im = plt.pcolormesh(ac, ac, Z_t, cmap=cmap, norm=norm, rasterized=True)
            ax.set_xscale("log")
            ax.set_yscale("log")
            plt.axis("scaled")

            t = format_seconds_as_years(tc[i_t])
            plt.title("$t = " + f"{t}" + "$")
            if subplot_idx + 3 >= len(times_of_interest):
                plt.xlabel("Particle Radius $a_j$ [m]")
            if subplot_idx % 3 == 0:
                plt.ylabel("Particle Radius $a_i$ [m]")

            ax.set_xscale("log")
            ax.set_yscale("log")
            plt.axis("scaled")
            plt.tight_layout()

        if symbol == "P_ij":
            label = "$P_{ij}^{sample}$"
        elif symbol == "S_ij":
            label = "$N_{ij}^{sample}$"
        else:
            raise Exception("")

        plt.subplots_adjust(top=0.9)
        # cax = fig.add_axes([0.93, 0.2, 0.013, 0.6])
        cax = fig.add_axes([0.1, 0.95, 0.8, 0.02])
        # plt.colorbar(cax=cax, orientation="vertical")
        cb = plt.colorbar(cax=cax, orientation="horizontal")
        cb.ax.set_title(label)

        plt.savefig(Path(path_to_figures, f"{title} {symbol}.pdf"))
        if should_show_plot:
            plt.show()
        plt.close()


if __name__ == "__main__":
    cfgs = []

    SAMPLING_DENSITIES = [0.5]
    N_m = 50

    cfgs, titles = [], []
    for rho_sample in SAMPLING_DENSITIES:
        for (enable_coagulation, enable_fragmentation) in [(True, True), (True, False)]: # (False, True)
            nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
            cfg = Config(
                mass_resolution=N_m,
                enable_collision_sampling=True,
                enable_coagulation=enable_coagulation,
                enable_fragmentation=enable_fragmentation,
                nr_of_samples=nr_of_samples,
            )
            cfgs.append(cfg)
            title = f"coag={enable_coagulation} frag={enable_fragmentation} {rho_sample=}"
            titles.append(title)

    for cfg, title in tqdm(zip(cfgs, titles)):
        main(cfg, title, should_show_plot=True)
