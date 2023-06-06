import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from utils.errors import handle_unknown_scale

SLIDER_COORDS = [0.925, 0.3, 0.04, 0.4]  # todo: move elsewhere ?
FIGSIZE = (6, 6)


class InteractiveSliderLinePlot:

    def __init__(
        self,
        cfg,
        x,
        ys,
        i_y=0,
        figsize=None, title=None,
        xlims=(None, None), ylims=(None, None),
        label=None, xlabel=None, ylabel=None,
    ):
        self.cfg = cfg
        self.x = x
        self.ys = ys
        self.i_y = i_y

        self.title = title
        self.xlims = xlims
        self.ylims = ylims
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.label = label

        if figsize is None:
            figsize = FIGSIZE
        self.fig = plt.figure(figsize=figsize)
        self.ax_plot = self.fig.add_subplot(111)
        self.ax_slider = self.fig.add_axes(SLIDER_COORDS)

        valmax = self.ys.shape[0] - 1
        self.slider = Slider(
            ax=self.ax_slider, label="", orientation="vertical",
            valmin=0, valmax=valmax, valstep=1, valinit=i_y
        )
        self.slider.on_changed(self.update)

    def draw(self):
        self.ax_plot.clear()
        self.fig.canvas.draw_idle()
        if self.cfg.mass_axis_scale == "lin":
            self.ax_plot.plot(self.x, self.ys[self.i_y], label=self.label)
        elif self.cfg.mass_axis_scale == "log":
            self.ax_plot.loglog(self.x, self.ys[self.i_y], label=self.label)
        else:
            handle_unknown_scale(self.cfg.mass_axis_scale)
        self.set_limits()
        self.set_texts()
        # plt.legend() # TODO Why not working?

    def update(self, i_y):
        self.i_y = i_y
        self.draw()

    def set_limits(self):
        if self.xlims[0] is not None:
            self.ax_plot.set_xlim(left=self.xlims[0])
        if self.xlims[1] is not None:
            self.ax_plot.set_xlim(right=self.xlims[1])
        if self.ylims[0] is not None:
            self.ax_plot.set_ylim(bottom=self.ylims[0])
        if self.ylims[1] is not None:
            self.ax_plot.set_ylim(top=self.ylims[1])

    def set_texts(self):
        if self.title is not None:
            self.ax_plot.set_title(self.title)
        if self.xlabel is not None:
            self.ax_plot.set_xlabel(self.xlabel)
        if self.ylabel is not None:
            self.ax_plot.set_ylabel(self.ylabel)
