from matplotlib import colors
from matplotlib.widgets import Slider
import matplotlib.pyplot as plt
import numpy as np

from visualization.kernel.config import FIGSIZE_1, FIGSIZE_2, FIGSIZE_3
from visualization.kernel.config import SLIDER_GEOMETRY, COLORBAR_GEOMETRY


class InteractiveKernelLayerPlot():

    def __init__(
        self, 
        kernels,
        k=25,
        figsize=None,
        symmetrize_kernels=False,
        kernel_subplot_titles=[],
        cmap_limits=None,
    ):
        # Pre-define `figsize` for 1, 2, or 3 kernels.
        if figsize is None:
            if len(kernels) == 1:
                figsize = FIGSIZE_1
            elif len(kernels) == 2:
                figsize = FIGSIZE_2
            elif len(kernels) == 3:
                figsize = FIGSIZE_3
            else:
                raise Exception("Please define & pass `figsize` argument.")
        # Define figure with axes for slider, kernel plots, & their colorbars.
        self.fig = plt.figure(figsize=figsize)
        # Define axis for slider.
        self.slider_ax = self.fig.add_axes(SLIDER_GEOMETRY)
        # Define axes for kernel plots.
        self.kernel_axes = [
            plt.subplot(1, len(kernels), i+1) for i, _ in enumerate(kernels)
        ]
        # Define axes (& their dimensions) for kernel colorbars.
        self.colorbar_axes = [
            self.fig.add_axes(COLORBAR_GEOMETRY) for _ in kernels
        ]
        # Define field for list of kernel & make sure they have same dimensions.
        assert_equal_kernel_shapes(kernels)
        self.kernels = kernels
        # Define colormap norm for colorbar.
        if cmap_limits is None:
            vmax = max([float(np.abs(K).max()) for K in self.kernels])
            vmin = -vmax
        else:
            vmin = cmap_limits[0]
            vmax = cmap_limits[1]
        self.cmap_norm = colors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
        # Define field for list of subplot titles.
        self.kernel_subplot_titles = kernel_subplot_titles
        # Define initial `k` value ("kernel layer"), controlled by slider.
        self.k = k
        # Setup the slider & initialize display of plot.
        self.symmetrize_kernels = symmetrize_kernels  # TODO Define a "kernel subplot" class?
        self.setup_slider()
        self.update(None)

    def setup_slider(self):
        # Get kernel dimension for max. value of slider.
        N_m = self.kernels[0].shape[0]  # <- See `assert_equal_kernel_shapes`.
        # Define slider & its `on_changed` method.
        self.slider = Slider(
            ax=self.slider_ax, label="$k$", orientation="horizontal",
            valmin=0, valmax=N_m-1, valstep=1, valinit=self.k,
        )
        self.slider.on_changed(self.update)

    def update(self, k):
        is_showing_for_first_time = k is None
        self.k = int(self.slider.val)
        for ax_idx in range(len(self.kernel_axes)):
            self.plot(ax_idx, self.k, is_showing_for_first_time)

    def plot(self, ax_idx, k, is_showing_for_first_time):
        # Clear axis of kernel plot. -> Better performance.
        ax = self.kernel_axes[ax_idx]
        ax.clear()
        plt.sca(ax)
        # Get kernel & kernel layer.
        K = self.kernels[ax_idx]
        K_k = K[k]
        if self.symmetrize_kernels:
            K_k = 0.5 * (K_k + K_k.T)
        # Draw plot.
        im = plt.pcolor(K_k, norm=self.cmap_norm)
        plt.axis("scaled")
        # Configure plot (e.g. labels).
        try:
            title = self.kernel_subplot_titles[ax_idx]
            plt.title(f"{title}" + " for $k=" + f"{k}" + "$")
        except IndexError:
            print(f"WARNING: No label found for plot on axis {ax_idx}.")
        plt.ylabel("$i$", rotation=0)
        plt.xlabel("$j$")
        # Draw colorbar, but only on first time shown. -> Better performance.
        plt.set_cmap("bwr")
        ax = self.colorbar_axes[ax_idx]
        if is_showing_for_first_time:
            self.fig.colorbar(im, cax=ax, orientation='horizontal')

    def show(self):
        plt.show()
        plt.close()  # Close plot. -> Better performance (saving RAM).


def assert_equal_kernel_shapes(kernels):
    if len(kernels) == 0:
        raise Exception("ERROR: You should give at least 1 kernel.")
    shape = kernels[0].shape
    for kernel in kernels[1:]:
        if kernel.shape != shape:
            raise Exception("ERROR: Kernel shapes do not match.")
    if len(shape) == 0:
        raise Exception("Error: Kernel shape has length 0.")
    dim = shape[0]
    for d in shape[1:]:
        if d != dim:
            raise Exception("ERROR: Kernel ist not cubic in shape.")
