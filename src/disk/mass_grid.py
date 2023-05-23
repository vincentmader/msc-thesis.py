from utils.axis import DiscreteAxis


class MassGrid(DiscreteAxis):
    def __init__(self, cfg):
        super().__init__(
            cfg.mass_min_value,
            cfg.mass_max_value,
            cfg.mass_resolution,
            cfg.mass_axis_scale,
        )
