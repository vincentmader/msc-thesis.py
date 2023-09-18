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
        """
        self.x_min = x_min
        self.x_max = x_max
        self.N = N
        self.scale = scale

        self.bin_boundaries = self._bin_boundaries()
        self.bin_centers = self._bin_centers()
        self.bin_widths = self._bin_widths()

    def _bin_boundaries(self) -> np.ndarray:
        x_min, x_max, N = self.x_min, self.x_max, self.N
        indices = np.arange(0, N + 1, 1)
        if self.scale == "lin":
            x = x_min + (x_max - x_min) * (indices / N)
            assert x[0] == x_min and x[N] == x_max
            return x
        if self.scale == "log":
            x = x_min * (x_max / x_min) ** (indices / N)
            assert x[0] == x_min and x[N] == x_max
            return x
        raise Exception(f"Axis scale '{self.scale}' unknown.")

    def _bin_centers(self) -> np.ndarray:
        xs = self.bin_boundaries
        if self.scale == "lin":
            return (xs[:-1] + xs[1:]) / 2
        if self.scale == "log":
            return np.sqrt(xs[:-1] * xs[1:])
        raise Exception(f"Axis scale '{self.scale}' unknown.")

    def _bin_widths(self) -> np.ndarray:
        bin_boundaries = self.bin_boundaries
        return bin_boundaries[1:] - bin_boundaries[:-1]

    def index_from_value(self, x) -> int:
        xs = self.bin_centers
        i = np.where((xs - x) <= 0)[0][-1]
        # NOTE: I believe this automatically prevents an `i >= N`.
        return i

    def value_from_index(self, i) -> np.float64:
        xs = self.bin_centers
        x = xs[i]
        return np.float64(x)
