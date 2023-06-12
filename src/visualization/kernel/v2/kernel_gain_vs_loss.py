import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


FIGSIZE = (12, 6)
SLIDER_RECT = [0.05, 0.5, 0.01, 0.5]
COLORBAR_GEOMETRY = [0.5, 0.95, 0.3, 0.01]
LOG_CMAP_VMIN = 1e-14


class KernelGainVsLossPlot:

    def __init__(
        self,
        kernel,
        enable_slider=True,
        slider_orientation="vertical",
        k=None,
    ):
        self.is_showing_for_the_first_time = True

        self.kernel = kernel
        self.mg = kernel.mg
        self.scale = kernel.mg.scale

        k = k if k else self.mg.N // 2

        self.fig = plt.figure(figsize=FIGSIZE)
        self.ax_gain = plt.subplot(1, 2, 1)
        self.ax_loss = plt.subplot(1, 2, 2)

        self.cmap_limits = cmap_limits(kernel)

        if enable_slider:
            self.setup_slider(k, slider_orientation)

        self.draw(k)

    def draw(self, k):
        # Draw gain & loss terms separately.
        p_gain = self.draw_gain(k)
        p_loss = self.draw_loss(k)
        # Setup colorbar.
        if self.is_showing_for_the_first_time:
            for p_id, p in enumerate([p_gain, p_loss]):
                rect = COLORBAR_GEOMETRY
                rect[0] += 0.5 * p_id - 0.5 * rect[2]
                self.setup_colorbar(p, rect, "horizontal")

    def draw_gain(self, k):
        K_gain = self.kernel.K_gain
        plt.sca(self.ax_gain)
        title = "K_{kij}^{gain}"
        return self.draw_kernel(K_gain, k, title)

    def draw_loss(self, k):
        K_loss = self.kernel.K_loss
        K_loss = -K_loss if self.scale == "log" else K_loss
        plt.sca(self.ax_loss)
        title = "K_{kij}^{loss}"
        title = f"-{title}" if self.scale == "log" else title
        return self.draw_kernel(K_loss, k, title)

    def draw_kernel(self, K, k, title):
        scale = self.kernel.mg.scale
        assert scale in ["lin", "log"], f"Invalid scale: {scale}"
        K_k = K[k]
        K_k = 0.5 * (K_k + K_k.T)
        K_k[K_k==0] = 1e-20 # TODO
        title = title if scale == "lin" else f"\log({title})"
        title = f"${title}$"
        plt.title(title)
        ax = plt.gca()
        ax.set_aspect('equal', adjustable='box')
        vmin, vmax = self.cmap_limits
        if scale == "log":
            plt.set_cmap("Reds")
            cmap_norm = colors.LogNorm(vmin=vmin, vmax=vmax)
        else:
            plt.set_cmap("bwr")
            cmap_norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
        return plt.pcolor(K_k, norm=cmap_norm)

    def clear(self):
        self.ax_gain.clear()
        self.ax_loss.clear()

    def show(self):
        plt.show()
        plt.close()

    def save(self, path):
        plt.savefig(path)

    def setup_colorbar(self, p, rect, orientation):
        rect = center_rect(rect, orientation)
        ax = self.fig.add_axes(rect)
        self.fig.colorbar(p, cax=ax, orientation=orientation)

    def setup_slider(self, k, orientation):
        from matplotlib.widgets import Slider
        rect = center_rect(SLIDER_RECT, orientation)
        ax = self.fig.add_axes(rect)
        slider = Slider(
            ax=ax, orientation=orientation, label="$k$", # todo: Make label customizable?
            valmin=0, valmax=self.mg.N-1, valstep=1, valinit=k,
        )
        slider.on_changed(self.slider_on_changed)
        self.slider = slider

    def slider_on_changed(self, slider_val):
        self.is_showing_for_the_first_time = False
        self.clear()
        self.draw(slider_val)


def center_rect(rect, orientation): 
    assert orientation in ["vertical", "horizontal"], f"Invalid orientation: {orientation}"
    if orientation == "vertical":
        return [rect[0] - rect[2], rect[1] - rect[3]/2, rect[2]*2, rect[3]]
    if orientation == "horizontal":
        return [rect[0] - rect[2]/2, rect[1] - rect[3], rect[2], rect[3]*2]

def cmap_limits(kernel) -> list[float]:
    scale = kernel.mg.scale

    K_gain = kernel.K_gain
    K_loss = kernel.K_loss

    if scale == "lin":
        vmin = min([K_gain.min(), K_loss.min()])
        vmax = max([K_gain.max(), K_loss.max()])
    elif scale == "log":
        K_loss = -K_loss # TODO needed?
        vmin = LOG_CMAP_VMIN
        vmax = max([K_gain.max(), K_loss.max()])
    else:
        raise Exception(f"Invalid scale: {scale}")

    if vmin == vmax:
        return [-1, 1]
    return [vmin, vmax]
