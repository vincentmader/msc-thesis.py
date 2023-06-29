from .axis import DiscreteAxis


class RadialGrid(DiscreteAxis):

    def __init__(self, cfg):
        super().__init__(
            cfg.radial_min_value,
            cfg.radial_max_value,
            cfg.radial_resolution,
            cfg.radial_axis_scale,
        )
