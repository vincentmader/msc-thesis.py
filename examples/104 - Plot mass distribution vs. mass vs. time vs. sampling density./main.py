import os, sys
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from models.solver import SolverV2
from constants import SECONDS_PER_YEAR

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "104")
path_to_figures  = Path(PATH_TO_FIGURES, "104")

# TS = list(range(155, 190, 1))
TS = list(range(150, 190, 2))


def integrate(
    cfg: Config,
    path: Path,
    to_save: list[str]
):
    solver = SolverV2(cfg)
    solver.run()
    
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

def main(should_show_plot=False):
    to_save = ["tc", "mb_k", "mc_k", "dm_k", "n_k", "N_k", "M_k", "dMdt_k"]
    for N_m in [50]:
        for enable_coagulation in [True, False]:
            for enable_fragmentation in [True, False]:
                if enable_coagulation == enable_fragmentation == False:
                    continue

                plt.figure(figsize=(10, 11))

                sampling_densities = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0]
                enable_collision_sampling = True
                for subplot_idx, rho_sample in enumerate(sampling_densities):
                    nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)

                    if  enable_fragmentation == True and enable_coagulation == False:
                        initial_mass_bin = 40
                    else:
                        initial_mass_bin = 0

                    cfg = Config(
                        mass_resolution=N_m,
                        enable_collision_sampling=True if rho_sample < 1.0 else False,
                        enable_coagulation=enable_coagulation,
                        enable_fragmentation=enable_fragmentation,
                        nr_of_samples=nr_of_samples,
                        initial_mass_bin=initial_mass_bin,
                    )

                    foo = lambda x: "enabled" if x else "disabled"
                    title  = "mass distribution with"
                    title += f" sampling {foo(enable_collision_sampling)},"
                    title += f" {rho_sample=}," # if enable_collision_sampling else ""
                    title += f" coag. {foo(enable_coagulation)},"
                    title += f" frag. {foo(enable_fragmentation)},"
                    title += f" {N_m=},"

                    path = Path(path_to_outfiles, title)
                    if not os.path.exists(path):
                        options = integrate(cfg, path, to_save)
                    else:
                        tc     = np.loadtxt(Path(path_to_outfiles, title, "tc.txt"))
                        mb_k   = np.loadtxt(Path(path_to_outfiles, title, "mb_k.txt"))
                        mc_k   = np.loadtxt(Path(path_to_outfiles, title, "mc_k.txt"))
                        dm_k   = np.loadtxt(Path(path_to_outfiles, title, "dm_k.txt"))
                        n_k    = np.loadtxt(Path(path_to_outfiles, title, "n_k.txt"))
                        N_k    = np.loadtxt(Path(path_to_outfiles, title, "N_k.txt"))
                        M_k    = np.loadtxt(Path(path_to_outfiles, title, "M_k.txt"))
                        dMdt_k = np.loadtxt(Path(path_to_outfiles, title, "dMdt_k.txt"))
                        options = [tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k]

                    ax = plt.subplot(4, 2, subplot_idx+1)
                    tc, mb_k, mc_k, dm_k, n_k, N_k, M_k, dMdt_k = options

                    colors = plt.cm.YlGnBu(np.linspace(0, 0.9, len(TS)))
                    for ii, i_t in enumerate(TS):
                        t = i_t # TODO
                        from functions.utils.dates import format_seconds_as_years
                        t = format_seconds_as_years(tc[i_t])
                        plt.loglog(mc_k, M_k[i_t], label=f"t={t}", color=colors[ii])
                        plt.title(r"$\rho_{sample}=$" + f"{rho_sample}")
                        # plt.ylim(1e-27, 1e-19)
                        plt.ylim(1e-17, 1e-7)
                    if subplot_idx % 2 != 0:
                        # plt.yticks([])
                        pass
                    else:
                        plt.ylabel(r"Dust Density $\rho_i^d$ [kg m$^{-3}$]")
                    if subplot_idx < len(sampling_densities) - 2:
                        ax.set_xticklabels([])
                    else:
                        plt.xlabel("Dust Particle Mass $m^c_i$ [kg]")
                    plt.xlim(mb_k[0], mb_k[-1])
                    plt.grid(True)

                # plt.tight_layout()

                cax = plt.axes([0.13, 0.95, 0.77, 0.02])
                cb = mpl.colorbar.ColorbarBase(cax, 
                    orientation="horizontal", 
                    cmap="YlGnBu",
                    norm=mpl.colors.LogNorm(
                        tc[TS[0]] / SECONDS_PER_YEAR, 
                        tc[TS[-1]] / SECONDS_PER_YEAR,
                    ),
                )
                cb.ax.set_title("Time $t$ [y]")

                os.makedirs(path_to_figures, exist_ok=True)
                name = f"3x2 rho_d vs. m, t, rho_sample, {N_m=}, coag={enable_coagulation}, frag={enable_fragmentation}.pdf"
                plt.savefig(Path(path_to_figures, name))
                if should_show_plot:
                    plt.show()
                plt.close()

if __name__ == "__main__":
    main(should_show_plot=False)
