import numpy as np
from utils.errors import handle_unknown_scale


class DiscreteAxis:
    def __init__(self, x_min, x_max, N_x, scale):
        self.x_min = x_min
        self.x_max = x_max
        self.N_x = N_x
        self.scale = scale

    def indices(self):
        return np.arange(0, self.N_x, 1)

    def grid_cell_boundaries(self):
        x_min, x_max, N_x = self.x_min, self.x_max, self.N_x
        xs = np.linspace(0, 1, N_x + 1)
        if self.scale == "lin":
            return x_min + (x_max - x_min) * xs
        if self.scale == "log":
            return x_min * (x_max / x_min) ** xs
        handle_unknown_scale(self.scale)

    def grid_cell_centers(self):
        x = self.grid_cell_boundaries()
        if self.scale == "lin":
            return (x[:-1] + x[1:]) / 2
        if self.scale == "log":
            return np.sqrt(x[:-1] * x[1:])
        handle_unknown_scale(self.scale)

    def grid_cell_width(self) -> float:
        x_min, x_max, N_x = self.x_min, self.x_max, self.N_x
        if self.scale == "lin":
            return (x_max - x_min) / N_x
        if self.scale == "log":
            return (x_max / x_min)**(1 / N_x)
        handle_unknown_scale(self.scale)

    def index_from_value(self, x):
        x_min, dx = self.x_min, self.grid_cell_width()
        if self.scale == "lin":
            res = (x - x_min) / dx
            return res.astype(int)
        if self.scale == "log":
            res = np.log(x / x_min) / np.log(dx)
            return res.astype(int)
        handle_unkown_scale(self.scale)

    def value_from_index(self, i):
        x_min, dx = self.x_min, self.grid_cell_width()
        if self.scale == "lin":
            return x_min + dx * i
        if self.scale == "log":
            return x_min * dx**i
        handle_unkown_scale(self.scale)
