from typing import Optional

import matplotlib.pyplot as plt


class BasePlot:
    __slots__ = ["fig"]

    def __init__(
        self,
        figsize: Optional[tuple[int, int]] = None,
    ):
        self.fig = plt.figure(figsize=figsize)

    def render(
        self,
        show_plot=True,
        save_plot=False,
        path_to_outfile=None,
    ):
        self.draw()
        if show_plot:
            plt.show()
        if save_plot:
            assert path_to_outfile is not None, "No savefile path specified."
            plt.savefig(path_to_outfile)
        plt.close()

    def draw(self):
        raise Exception("Not implemented.")
