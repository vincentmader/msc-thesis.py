import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_OUTFILES
from models.solver import SolverV2


def integrate(
    cfg: Config,
    title: str,
    to_save: list[str]
):
    solver = SolverV2(cfg)
    solver.run()
    
    path = Path(PATH_TO_OUTFILES, "data", title)
    solver.save(path, to_save)
    return [
        solver.tg.bin_centers,
        solver.mg.bin_boundaries,
        solver.mg.bin_centers,
        solver.mg.bin_widths,
        solver.n_k_vs_t,
        solver.N_k_vs_t,
        solver.M_k_vs_t,
        solver.dMdt_vs_t,
    ]

def has_integrated(title: str):
    if os.path.exists(Path(PATH_TO_OUTFILES, "data", title)):
        return True
    return False

def plot_1(
    title,
    options,
):
    tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k = options

    TS = list(range(155, 190, 1))
    # TS = list(range(0, 120, 1))

    M_k = N_k * mc_k * dm_k
    plt.figure(figsize=(11, 4))
    plt.ylim(1e-32, 1e-16)
    plt.title(title)
    # plt.title(f"${N_m=}$, " + r"$\rho_{sample}=$" + f"{rho_sample}")
    # path = Path(PATH_TO_OUTFILES, title, f"{N_m=}", f"{rho_sample=}", "M_k.pdf")
    # plot(M_k, times_of_interest, mb_k, mc_k, tc, path=path)

    plt.xlim(mb_k[0], mb_k[-1]*20)

    # TS = list(range(0, 155, 1))
    # colors = plt.cm.jet(np.linspace(0,1,len(TS)))
    colors = plt.cm.Reds(np.linspace(0,1,len(TS)))
    for ii, i_t in enumerate(TS):
        t = i_t # TODO
        from functions.utils.dates import format_seconds_as_years
        t = format_seconds_as_years(tc[i_t])
        plt.loglog(mc_k, M_k[i_t], label=f"t={t}", color=colors[ii])

    plt.xlabel("mass $m_i$")
    plt.ylabel(r"$N_i \cdot m_i \cdot \Delta m_i$")
    plt.legend(loc="upper right")

    # if path:
    #     plt.savefig(path)
    plt.show()
    plt.close()


def plot_2(title, options):
    tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k = options

    TS = list(range(155, 190, 1))
    colors = plt.cm.Reds(np.linspace(0, 1, len(TS)))
    for ii, i_t in enumerate(TS):
        plt.loglog(mc_k, dMdt_k[i_t], color=colors[ii])

    plt.show()
    plt.close()

def plot_3(title, options):
    tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k = options

    M = np.sum(M_k, axis=1)
    err = (M - M[0]) / M[0] * 100
    x, y = tc, err

    plt.figure(figsize=(11, 4))
    plt.loglog(x, y, "r")
    plt.loglog(x, -y, "b")
    plt.title(title)

    plt.ylim(0, 1e-10)
    plt.show()
    plt.close()


def v1():
    
    to_save = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k"]

    for N_m in [25, 50]:
        for enable_coagulation in [True]:
            for enable_fragmentation in [True, False]:
                # for enable_collision_sampling in [True, False]:
                enable_collision_sampling = True
                for rho_sample in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:

                    initial_mass_bin = 49 if (enable_fragmentation and not enable_coagulation) else 0

                    nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)

                    cfg = Config(
                        mass_resolution=N_m,
                        # enable_collision_sampling=enable_collision_sampling,
                        enable_collision_sampling=True if rho_sample < 1.0 else False,
                        enable_coagulation=enable_coagulation,
                        enable_fragmentation=enable_fragmentation,
                        initial_mass_bin=initial_mass_bin,
                        nr_of_samples=nr_of_samples,
                    )

                    foo = lambda x: "enabled" if x else "disabled"
                    title  = "mass distribution with"
                    title += f" sampling {foo(enable_collision_sampling)},"
                    title += f" {rho_sample=}," # if enable_collision_sampling else ""
                    title += f" coag. {foo(enable_coagulation)},"
                    title += f" frag. {foo(enable_fragmentation)},"
                    title += f" {N_m=},"

                    if not has_integrated(title):
                        options = integrate(cfg, title, to_save)
                    else:
                        tc     = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "tc.txt"))
                        mb_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "mb_k.txt"))
                        mc_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "mc_k.txt"))
                        dm_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "dm_k.txt"))
                        n_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "n_k.txt"))
                        N_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "N_k.txt"))
                        M_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "M_k.txt"))
                        dMdt_k = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "dMdt_k.txt"))
                        options = [tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k]

                    # plot_1(title, options)  # N*m*dm vs. t
                    # plot_2(title, options)  # dN/dt*m*dm vs. t
                    plot_3(title, options)  # err_M vs. t


def v2():
    to_save = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k"]
    for N_m in [25, 50]:
        for enable_coagulation in [True]:
            for enable_fragmentation in [True, False]:

                plt.figure(figsize=(7, 9))

                sampling_densities = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
                enable_collision_sampling = True
                for subplot_idx, rho_sample in enumerate(sampling_densities):
                    nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
                    cfg = Config(
                        mass_resolution=N_m,
                        enable_collision_sampling=True if rho_sample < 1.0 else False,
                        enable_coagulation=enable_coagulation,
                        enable_fragmentation=enable_fragmentation,
                        nr_of_samples=nr_of_samples,
                    )

                    foo = lambda x: "enabled" if x else "disabled"
                    title  = "mass distribution with"
                    title += f" sampling {foo(enable_collision_sampling)},"
                    title += f" {rho_sample=}," # if enable_collision_sampling else ""
                    title += f" coag. {foo(enable_coagulation)},"
                    title += f" frag. {foo(enable_fragmentation)},"
                    title += f" {N_m=},"

                    if not has_integrated(title):
                        options = integrate(cfg, title, to_save)
                    else:
                        tc     = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "tc.txt"))
                        mb_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "mb_k.txt"))
                        mc_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "mc_k.txt"))
                        dm_k   = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "dm_k.txt"))
                        n_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "n_k.txt"))
                        N_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "N_k.txt"))
                        M_k    = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "M_k.txt"))
                        dMdt_k = np.loadtxt(Path(PATH_TO_OUTFILES, "data", title, "dMdt_k.txt"))
                        options = [tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k]

                    plt.subplot(4, 2, subplot_idx+1)
                    tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k = options
                    TS = list(range(155, 190, 1))
                    M_k = N_k * mc_k * dm_k

                    colors = plt.cm.Blues(np.linspace(0, 1, len(TS)))
                    for ii, i_t in enumerate(TS):
                        t = i_t # TODO
                        from functions.utils.dates import format_seconds_as_years
                        t = format_seconds_as_years(tc[i_t])
                        plt.loglog(mc_k, M_k[i_t], label=f"t={t}", color=colors[ii])
                        plt.title(r"$\rho_{sample}=$" + f"{rho_sample}")
                        plt.ylim(1e-33, 1e-18)
                    if subplot_idx % 2 != 0:
                        plt.yticks([])
                    else:
                        plt.ylabel(r"$n_k \cdot m_k \cdot \Delta m_k$")
                    if subplot_idx < len(sampling_densities) - 2:
                        plt.xticks([])
                    else:
                        plt.xlabel("time")

                plt.tight_layout()
                path = Path(PATH_TO_OUTFILES, "data")
                os.makedirs(path, exist_ok=True)
                plt.savefig(Path(path, f"3x2 rho_d vs. m, t, rho_sample, {N_m=}, {enable_fragmentation=}.pdf"))
                plt.show()
                plt.close()


v1()
# v2()
