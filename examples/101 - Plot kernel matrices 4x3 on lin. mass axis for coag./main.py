import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import PATH_TO_FIGURES
from models.kernel import Kernel

path_to_figures = Path(PATH_TO_FIGURES, "101")
os.makedirs(path_to_figures, exist_ok=True)


N_m = 100
ks = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96]

vmin = -0.5
vmax = +0.5


def create_4x3_plot(
    cfg, 
    should_show_plot=False, 
    should_symmetrize_kernel_layers=False
):

    kernel = Kernel(cfg)
    mc = kernel.mg.bin_centers
    K = kernel.K / kernel.R_coll  # <- Normalize kernel for linear plot.
    
    fig = plt.figure(figsize=(10, 11))
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)

    for subplot_idx, k in enumerate(ks):

        K_k = K[int(k)]
        if should_symmetrize_kernel_layers:
            K_k = (K_k + K_k.T) / 2

        ax = plt.subplot(4, 3, subplot_idx + 1)
        im = plt.pcolormesh(mc, mc, K_k, cmap="bwr", norm=norm, rasterized=True)

        plt.title("$"+ f"{k=}" + "$")
        if subplot_idx + 3 >= len(ks):
            plt.xlabel("particle mass $m_j$ [kg]")
        if subplot_idx % 3 == 0:
            plt.ylabel("particle mass $m_i$ [kg]")

        plt.axis("scaled")
        plt.tight_layout()

    cax = fig.add_axes([0.93, 0.2, 0.013, 0.6])
    plt.colorbar(cax=cax, orientation="vertical")

    name = f"Kkij vs k, coag={cfg.enable_coagulation}, frag={cfg.enable_fragmentation}.pdf"
    plt.savefig(Path(path_to_figures, name))
    if should_show_plot:
        plt.show()
    plt.close()


def main(
    should_show_plot=False,
    should_symmetrize_kernel_layers=False,
):
    for (enable_coagulation, enable_fragmentation) in tqdm([(True, True), (True, False), (False, True)]):
    
        cfg = Config(
            mass_min_value=1,
            mass_max_value=100,
            mass_resolution=100,
            mass_axis_scale="lin",
            enable_coagulation=enable_coagulation,
            enable_fragmentation=enable_fragmentation,
        )
    
        create_4x3_plot(
            cfg, 
            should_show_plot=should_show_plot,
            should_symmetrize_kernel_layers=should_symmetrize_kernel_layers,
        )

if __name__ == "__main__":
    main(
        should_show_plot=False,
        should_symmetrize_kernel_layers=True,
    )
