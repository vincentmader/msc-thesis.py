from typing import Optional

import matplotlib.pyplot as plt


class BasePlot:
    __slots__ = ["fig", "xlims", "ylims"]

    def __init__(
        self,
        figsize: Optional[tuple[int, int]] = None,
        xlims: Optional[tuple[int, int]] = None,
        ylims: Optional[tuple[int, int]] = None,
    ):
        self.fig = plt.figure(figsize=figsize)
        self.xlims = xlims
        self.ylims = ylims

    def render(
        self,
        show_plot=True,
        save_plot=False,
        path_to_outfile=None,
    ):
        if self.xlims is not None:
            plt.xlim(self.xlims[0], self.xlims[1])
        if self.ylims is not None:
            plt.ylim(self.ylims[0], self.ylims[1])
        self.draw()
        if save_plot:
            assert path_to_outfile is not None, "No savefile path specified."
            plt.savefig(path_to_outfile)
        if show_plot:
            plt.show()
        plt.close()

    def draw(self):
        raise Exception("Not implemented.")
