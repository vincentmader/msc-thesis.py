from matplotlib.widgets import Slider

from .base_plot import BasePlot
from .gridspec_subplot import GridspecSubplot


class GridspecPlot(BasePlot):
    __slots__ = ["gs", "axes", "subplots", "slider"]

    def __init__(
        self,
        subplots: list[GridspecSubplot],
        add_slider=False,
    ):
        figsize = self._choose_figsize(len(subplots))
        super().__init__(figsize=figsize)
        self.subplots = subplots
        self._setup_axes(subplots, add_slider)

    def draw(self):
        self.fig.canvas.draw_idle()
        for idx, subplot in enumerate(self.subplots):
            axes = self.axes[2*idx: 2*idx+2]
            subplot.draw(axes)

    def _choose_figsize(
        self, 
        nr_of_subplots: int,
    ):
        figsizes = {
            1: (5, 6),
            2: (10, 6),
            3: (15, 6),
        }
        return figsizes[nr_of_subplots]

    def _setup_axes(
        self, 
        subplots: list[GridspecSubplot], 
        add_slider: bool,
    ):
        nr_of_rows = 2
        nr_of_cols = len(subplots) + 1 if add_slider else len(subplots)
        gs_dimensions = (nr_of_rows, nr_of_cols)

        gs_row_sizes = [1, 20]
        gs_col_sizes = [20] * len(subplots)
        gs_col_sizes += [1] if add_slider else []

        gs_locations = []
        for idx, subplot in enumerate(subplots):
            gs_location = (0, 1, idx, idx+1)
            gs_locations.append(gs_location)
            gs_location = (1, 2, idx, idx+1)
            gs_locations.append(gs_location)

        self.gs = self.fig.add_gridspec(
            nrows=gs_dimensions[0], 
            ncols=gs_dimensions[1],
            width_ratios=gs_col_sizes,
            height_ratios=gs_row_sizes,
            hspace=0.0,
            wspace=0.3,
        )

        self.axes = []
        for gs_location in gs_locations:
            row_i, row_f, col_i, col_f = gs_location
            foo = self.gs[row_i:row_f, col_i:col_f]
            ax = self.fig.add_subplot(foo)
            self.axes.append(ax)

        if add_slider:
            self._setup_slider()

    def _setup_slider(self):
        k_max = max([s.z.shape[2] - 1 for s in self.subplots]) # TODO
        k_init = k_max // 2 + 1

        ax = self.fig.add_subplot(self.gs[-1])
        self.slider = Slider(
            ax=ax, 
            valmin=0, 
            valmax=k_max, 
            valstep=1, 
            valinit=k_init,
            label="$k$", 
            orientation="vertical",
        )  # TODO
        self.slider.on_changed(self.update)
        self.axes.append(ax)

    def update(
        self, 
        k: int,
    ):
        for subplot in self.subplots:
            subplot.k = k
            subplot.update(k)
