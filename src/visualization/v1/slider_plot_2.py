import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


SLIDER_COORDS = [0.925, 0.3, 0.04, 0.4]
FIGSIZE = (10, 10)


class InteractiveSliderLinePlot:

    def __init__(
        self, cfg, x, ys_1, ys_2, i_y=0,
        figsize=None, title=None, label=None,
        xlims_1=(None, None), ylims_1=(None, None),
        xlims_2=(None, None), ylims_2=(None, None),
        xlabel_1=None, ylabel_1=None,
        xlabel_2=None, ylabel_2=None,
    ):
        # Save plot configuration class object fields.
        self.cfg = cfg
        self.x = x
        self.ys_1 = ys_1
        self.ys_2 = ys_2
        self.i_y = i_y
        self.title = title
        self.xlims_1 = xlims_1
        self.ylims_1 = ylims_1
        self.xlims_2 = xlims_2
        self.ylims_2 = ylims_2
        self.xlabel_1 = xlabel_1
        self.ylabel_1 = ylabel_1
        self.xlabel_2 = xlabel_2
        self.ylabel_2 = ylabel_2
        self.label = label

        # Setup figure.
        figsize = figsize if figsize is not None else FIGSIZE
        self.fig = plt.figure(figsize=figsize)

        # Setup two axes.
        self.ax_1 = self.fig.add_subplot(211)
        self.ax_slider = self.fig.add_axes(SLIDER_COORDS)
        self.ax_2 = self.fig.add_subplot(212)

        # Setup slider.
        valmax = self.ys_1.shape[0] - 1
        self.slider = Slider(
            ax=self.ax_slider, label="", orientation="vertical",
            valmin=0, valmax=valmax, valstep=1, valinit=i_y
        )
        self.slider.on_changed(self.update)

    def draw(self):
        self.fig.canvas.draw_idle()

        x = self.x
        y_1 = self.ys_1[self.i_y]
        y_2 = self.ys_2[self.i_y]

        x = x[:-1]
        y_1 = y_1[:-1]
        y_2 = y_2[:-1]

        # dt = self.x[1:] - self.x[:-1]   # NOTE Use this to plot dn/dt
        # dt = np.append(dt, dt[-1])      #      instead of dn/dt*DT
        # y_2 = y_2 / dt                      #

        # Plot particle mass distribution.
        if self.cfg.mass_axis_scale == "lin":
            self.ax_1.plot(x, y_1, label=self.label)
        elif self.cfg.mass_axis_scale == "log":
            self.ax_1.loglog(x, y_1, label=self.label)
        else:
            raise Exception(
                f"Axis scale '{self.cfg.mass_axis_scale}' unknown.")

        # Plot derivative of particle mass distribution.
        if self.cfg.mass_axis_scale == "lin":
            self.ax_2.plot(x, y_2)
        elif self.cfg.mass_axis_scale == "log":
            self.ax_2.loglog(x, y_2, 'r', label="pos.")
            self.ax_2.loglog(x, -y_2, 'b', label="neg.")
            self.ax_2.legend()
        else:
            raise Exception(
                f"Axis scale '{self.cfg.mass_axis_scale}' unknown.")

        self.set_labels()
        self.set_xlims()
        self.set_ylims()

    def update(self, i_y):
        self.i_y = i_y
        self.clear()
        self.draw()

    def clear(self):
        self.ax_1.clear()
        self.ax_2.clear()

    def set_xlims(self):
        self.ax_1.set_xlim(left=self.xlims_1[0])
        self.ax_1.set_xlim(right=self.xlims_1[1])
        self.ax_2.set_xlim(left=self.xlims_2[0])
        self.ax_2.set_xlim(right=self.xlims_2[1])

    def set_ylims(self):
        self.ax_1.set_ylim(bottom=self.ylims_1[0])
        self.ax_1.set_ylim(top=self.ylims_1[1])
        self.ax_2.set_ylim(bottom=self.ylims_2[0])
        self.ax_2.set_ylim(top=self.ylims_2[1])

    def set_labels(self):
        self.ax_1.set_title(self.title)
        self.ax_1.set_xlabel(self.xlabel_1)
        self.ax_1.set_ylabel(self.ylabel_1)
        self.ax_2.set_xlabel(self.xlabel_2)  # TODO Fix xticks.
        # TODO Draw real derivative.
        self.ax_2.set_ylabel(self.ylabel_2)
