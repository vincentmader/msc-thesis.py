from utils.axis import DiscreteAxis


class TimeGrid(DiscreteAxis):

    def __init__(self, cfg):
        super().__init__(
            cfg.time_min_value,
            cfg.time_max_value,
            cfg.time_resolution,
            cfg.time_axis_scale,
        )
