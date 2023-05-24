import numpy as np
from utils.errors import handle_unknown_scale


class DiscreteAxis:

    def __init__(self, x_min, x_max, N_x, scale):
        """
        Create a new `DiscreteAxis` object.

        Args:
            x_min:  This is the lower boundary of the discretized axis.
            x_max:  This is the upper boundary of the discretized axis.
            N_x:    This is the number of grid-points (bins) in the discretized axis.
            scale:  This should be either "lin" (linear) or "log" (logarithmic).

        Returns:
            None
        """
        self.x_min = x_min
        self.x_max = x_max
        self.N_x = N_x
        self.scale = scale

    def indices(self):
        return np.arange(0, self.N_x, 1)

    def grid_cell_boundaries(self) -> np.ndarray:
        x_min, x_max, N_x = self.x_min, self.x_max, self.N_x
        indices = np.linspace(0, 1, N_x + 1)
        if self.scale == "lin":
            return x_min + (x_max - x_min) * indices
        if self.scale == "log":
            return x_min * (x_max / x_min) ** indices
        handle_unknown_scale(self.scale)

    def grid_cell_centers(self) -> np.ndarray:
        xs = self.grid_cell_boundaries()
        if self.scale == "lin":
            return (xs[:-1] + xs[1:]) / 2
        if self.scale == "log":
            return np.sqrt(xs[:-1] * xs[1:])
        handle_unknown_scale(self.scale)

    def grid_cell_widths(self) -> np.ndarray:
        grid_cell_boundaries = self.grid_cell_boundaries()
        return grid_cell_boundaries[1:] - grid_cell_boundaries[:-1]

    def index_from_value(self, x) -> int:
        x_min, x_max, N_x = self.x_min, self.x_max, self.N_x
        if self.scale == "lin":
            d_x = (x_max - x_min) / N_x
            res = (x - x_min) / d_x
            return res.astype(int)
        if self.scale == "log":
            q_x = (x_max / x_min)**(1 / N_x)
            res = np.log(x / x_min) / np.log(q_x)
            return res.astype(int)
        handle_unkown_scale(self.scale)

    def value_from_index(self, i) -> float:
        x_min, x_max, N_x = self.x_min, self.x_max, self.N_x
        if self.scale == "lin":
            d_x = (x_max - x_min) / N_x
            return x_min + d_x * i
        if self.scale == "log":
            q_x = (x_max / x_min)**(1 / N_x)
            return x_min * q_x**i
        handle_unkown_scale(self.scale)
