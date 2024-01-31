import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from config import PATH_TO_OUTFILES
    from models.solver import SolverV2
except ModuleNotFoundError as e:
    raise e

""" Plot n*m*dm
    for different mass resolutions `N_m`
    and sampling densities `rho_sample`.
"""
MASS_GRID_RESOLUTIONS   = [50]
SAMPLING_DENSITIES      = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
TITLE                   = "2024-01-24"

def integrate():
    for N_m in MASS_GRID_RESOLUTIONS:
        for rho_sample in SAMPLING_DENSITIES:

            nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
            enable_collision_sampling = True if rho_sample < 1.0 else False

            cfg = Config(
                mass_resolution=N_m,
                nr_of_samples=nr_of_samples,
                enable_collision_sampling=enable_collision_sampling,
            )
            path = Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}")

            solver = SolverV2(cfg)
            solver.run()
            solver.save(path, ["N_k", "mb_k", "mc_k", "dm_k", "tc"])


def plot(y, times_of_interest, mb_k, mc_k, tc, path=None):
    plt.xlim(mb_k[0], mb_k[-1]*20)

    # TS = list(range(0, 155, 1))
    TS = times_of_interest
    colors = plt.cm.jet(np.linspace(0,1,len(TS)))
    colors = plt.cm.Blues(np.linspace(0,1,len(TS)))
    for ii, i_t in enumerate(TS):
        t = i_t # TODO
        from functions.utils.dates import format_seconds_as_years
        t = format_seconds_as_years(tc[i_t])
        plt.loglog(mc_k, y[i_t], label=f"t={t}", color=colors[ii])

    plt.xlabel("mass $m_i$")
    plt.ylabel(r"$N_i \cdot m_i \cdot \Delta m_i$")
    plt.legend(loc="upper right")

    if path:
        plt.savefig(path)
    plt.show()
    plt.close()


def f1():
    for N_m in MASS_GRID_RESOLUTIONS:
        for rho_sample in SAMPLING_DENSITIES:

            N_k  = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "N_k.txt"))
            mb_k = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "mb_k.txt"))
            mc_k = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "mc_k.txt"))
            dm_k = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "dm_k.txt"))
            tc   = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "tc.txt"))

            times_of_interest = list(range(155, 190, 1))

            M_k = N_k * mc_k * dm_k
            plt.figure(figsize=(15, 8))
            plt.ylim(1e-32, 1e-16)
            plt.title(f"${N_m=}$, " + r"$\rho_{sample}=$" + f"{rho_sample}")
            path = Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "M_k.pdf")
            plot(M_k, times_of_interest, mb_k, mc_k, tc, path=path)

            # N_k_complete = np.loadtxt(Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", "rho_sample=1.0", "N_k.txt"))
            # M_k_complete = N_k_complete * mc_k * dm_k
            # path = Path(PATH_TO_OUTFILES, TITLE, f"{N_m=}", f"{rho_sample=}", "M_k err rel to complete.pdf")
            # plt.figure(figsize=(15, 8))
            # plt.title(f"${N_m=}$, " + r"$\rho_{sample}=$" + f"{rho_sample}")
            # plot((M_k - M_k_complete) / M_k_complete, times_of_interest, mb_k, mc_k, tc, path=path)


def f2():
    TITLE = "2024-01-25_2"
    for enable_coagulation in [True, False]:
        for enable_fragmentation in [True, False]:
            for initial_mass_bin in [0, 40]:
                if enable_coagulation and initial_mass_bin != 0: 
                    continue # skip
                if not enable_coagulation and not enable_fragmentation:
                    continue
                cfg = Config(
                    enable_collision_sampling=False,
                    enable_coagulation=enable_coagulation,
                    enable_fragmentation=enable_fragmentation,
                    initial_mass_bin=initial_mass_bin
                )
                solver = SolverV2(cfg)
                solver.run()

                N_k = solver.N_k_vs_t
                mb_k = solver.mg.bin_boundaries
                mc_k = solver.mg.bin_centers
                dm_k = solver.mg.bin_widths
                tc = solver.tg.bin_centers

                times_of_interest = range(200)

                M_k = N_k * mc_k * dm_k
                plt.figure(figsize=(15, 8))
                plt.title(f"{enable_coagulation=}, {enable_fragmentation=}")
                plt.ylim(1e-35, 1e-15)
                path = Path(PATH_TO_OUTFILES, TITLE, f"M_k_vs_t {enable_coagulation=}, {enable_fragmentation=}.pdf")
                plot(M_k, times_of_interest, mb_k, mc_k, tc, path=path)


if __name__ == "__main__":

    if False:
        integrate()
    # f1()

    f2()
