import numpy as np

from .collision_cross_section import collision_cross_section
from functions.physics.dust.relative_velocity import relative_velocity


def collision_rate(cfg, disk, disk_region):  # TODO Think re: inputs of this function
    mg = disk.mg

    if cfg.enable_physical_collisions:
        sigma = collision_cross_section(cfg, mg)
        dv = relative_velocity(cfg, disk, disk_region)
        return sigma * dv     # NOTE Simplification was made here.

    # In the most simple case, the rates are just set to 1.
    else:  
        return np.ones(shape=[mg.N] * 2)
