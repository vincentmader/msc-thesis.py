import os, sys
from pathlib import Path
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES
from models.kernel import Kernel


FIGSIZE = (10, 11)

# TODOs:
# - TODO Implement saving to disk.

N_m = 100
ks = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96]
# ks = [10, 20, 30, 40, 50, 60, 70, 80, 90]

vmin = 1e-30
vmax = 1e-7


def create_4x3_plot(
    cfg, mc, ac, K, cmap, component, 
    should_show_plot=False, 
    should_symmetrize_kernel_layers=False
):

    fig = plt.figure(figsize=FIGSIZE)
    norm = LogNorm(vmin=vmin, vmax=vmax)

    for subplot_idx, k in enumerate(ks):

        K_k = K[int(k)]
        if should_symmetrize_kernel_layers:
            K_k = (K_k + K_k.T) / 2

        ax = plt.subplot(4, 3, subplot_idx + 1)
        im = plt.pcolormesh(
            ac, ac, K_k,
            cmap=cmap, norm=norm,
            rasterized=True,
        )

        plt.title("$"+ f"{k=}" + "$")
        if subplot_idx + 3 >= len(ks):
            plt.xlabel("Particle Radius $a_j$ [m]")
        if subplot_idx % 3 == 0:
            plt.ylabel("Particle Radius $a_i$ [m]")

        ax.set_xscale("log")
        ax.set_yscale("log")
        plt.axis("scaled")
        plt.tight_layout()

    label = "$K_{kij}^{" + component + "}$ [m$^3$ s$^{-1}$]"
    plt.subplots_adjust(top=0.9)
    # cax = fig.add_axes([0.93, 0.2, 0.013, 0.6])
    cax = fig.add_axes([0.1, 0.95, 0.8, 0.02])
    cb = plt.colorbar(cax=cax, orientation="horizontal")
    cb.ax.set_title(label)

    name = f"Kkij_{component} vs k, coag={cfg.enable_coagulation}, frag={cfg.enable_fragmentation}.pdf"
    path = Path(PATH_TO_FIGURES, "102")
    os.makedirs(path, exist_ok=True)
    plt.savefig(Path(path, name))
    if should_show_plot:
        plt.show()
    plt.close()


def main(
    should_show_plot=False,
    should_symmetrize_kernel_layers=False,
):
    cfgs = []
    
    enable_collision_sampling = False
    for (enable_coagulation, enable_fragmentation) in [(True, True), (True, False), (False, True)]:
        cfg = Config(
            mass_resolution=N_m,
            enable_collision_sampling=enable_collision_sampling,
            enable_coagulation=enable_coagulation, 
            enable_fragmentation=enable_fragmentation,
        )
        cfgs.append(cfg)
    
    for cfg in tqdm(cfgs):
        kernel = Kernel(cfg)
        mc = kernel.mg.bin_centers
        ac = kernel.mg.particle_radii

        for component in ["gain", "loss"]:

            if component == "gain":
                K = kernel.K_gain 
                cmap = "Reds"
            elif component == "loss":
                K = -kernel.K_loss 
                cmap = "Blues"
            else:
                raise Exception()

            create_4x3_plot(
                cfg, mc, ac, K, cmap, component,
                should_show_plot=should_show_plot,
                should_symmetrize_kernel_layers=should_symmetrize_kernel_layers,
            )

if __name__ == "__main__":
    main(                
        should_show_plot=False,
        should_symmetrize_kernel_layers=True,
    )
