import os
import sys
from pprint import pprint
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.widgets import Slider
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from kernel import Kernel
    from disk import MassGrid, Disk, RadialGrid, DiskRegion
    from utils.plotting import plt_show_then_close
    from config import Config
except ModuleNotFoundError as e:
    raise e

VMIN, VCEN, VMAX = -1, 0, 1
CMAP_NORM = colors.TwoSlopeNorm(vmin=VMIN, vcenter=VCEN, vmax=VMAX)
FIGSIZE = (6, 5)
SLIDER_POSITION = [0.05, 0.3, 0.02, 0.4]

cfg = Config(
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_physical_cross_sections=False,
    enable_physical_relative_velocities=[],
)
pprint(cfg.__dict__)

kernel = Kernel(cfg)
rg = RadialGrid(cfg)
mg = MassGrid(cfg)
disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
K = kernel.K(disk, disk_region)


def plot(k):
    plt.subplot(111)
    plt.title("$K_{kij}$ for $k=" + f"{k}" + "$")
    plt.set_cmap("bwr")
    K_k = K[k]
    K_k = 0.5 * (K_k + K_k.T)
    plt.pcolor(K_k, norm=CMAP_NORM)
    plt.colorbar()
    plt.axis("scaled")


def update(_):
    global k
    k = int(slider.val)
    plot(k)


if __name__ == "__main__":
    k = 25

    fig = plt.figure(figsize=FIGSIZE)
    ax_slider = fig.add_axes(SLIDER_POSITION)
    slider = Slider(
        ax=ax_slider,
        label="$k$",
        valmin=0,
        valmax=mg.N_x - 1,
        valstep=1,
        valinit=k,
        orientation="vertical"
    )
    slider.on_changed(update)
    fig.canvas.draw_idle()
    update(None)
    plt_show_then_close()
