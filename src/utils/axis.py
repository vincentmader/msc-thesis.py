import numpy as np


class DiscreteAxis:

    def __init__(self, x_min, x_max, N, scale):
        """
        Create a new `DiscreteAxis` object.

        Args:
            x_min:  This is the lower boundary of the discretized axis.
            x_max:  This is the upper boundary of the discretized axis.
            N:      This is the number of grid-points (bins) in the discretized axis.
            scale:  This should be either "lin" (linear) or "log" (logarithmic).

        Returns:
            None
        """
        self.x_min = x_min
        self.x_max = x_max
        self.N = N
        self.scale = scale

    def grid_cell_boundaries(self) -> np.ndarray:
        x_min, x_max, N = self.x_min, self.x_max, self.N
        indices = np.linspace(0, 1, N + 1)  # TODO Rename: These are not "indices".
        if self.scale == "lin":
            x = x_min + (x_max - x_min) * indices
            assert x[0] == x_min and x[-1] == x_max
            return x
        if self.scale == "log":
            x = x_min * (x_max / x_min) ** indices
            assert x[0] == x_min and x[-1] == x_max
            return x
        raise Exception(f"Axis scale '{self.scale}' unknown.")

    def grid_cell_centers(self) -> np.ndarray:
        xs = self.grid_cell_boundaries()
        if self.scale == "lin":
            return (xs[:-1] + xs[1:]) / 2
        if self.scale == "log":
            return np.sqrt(xs[:-1] * xs[1:])
        raise Exception(f"Axis scale '{self.scale}' unknown.")

    def grid_cell_widths(self) -> np.ndarray:
        grid_cell_boundaries = self.grid_cell_boundaries()
        return grid_cell_boundaries[1:] - grid_cell_boundaries[:-1]

    def index_from_value(self, x) -> int:
        x_min, x_max, N = self.x_min, self.x_max, self.N
        if self.scale == "lin":
            d_x = (x_max - x_min) / N
            res = (x - x_min) / d_x
            return res.astype(int)
        if self.scale == "log":
            q_x = (x_max / x_min)**(1 / N)
            res = np.log(x / x_min) / np.log(q_x)
            return res.astype(int)
        raise Exception(f"Axis scale '{self.scale}' unknown.")

    def value_from_index(self, i) -> np.float64:
        x_min, x_max, N = self.x_min, self.x_max, self.N
        if self.scale == "lin":
            d_x = (x_max - x_min) / N
            res = x_min + d_x * i
            return np.float64(res)
        if self.scale == "log":
            q_x = (x_max / x_min)**(1 / N)
            res = x_min * q_x**i
            return np.float64(res)
        raise Exception(f"Axis scale '{self.scale}' unknown.")
