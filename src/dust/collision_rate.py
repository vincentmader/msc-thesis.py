from .collision_cross_section import collision_cross_section
from .relative_velocity import relative_velocity


def collision_rate(cfg, disk, disk_region):  # TODO Think re: inputs of this function
    mg = disk.mass_axis

    # Calculate collision cross section & relative particle velocities.
    sigma = collision_cross_section(cfg, mg)
    dv = relative_velocity(cfg, disk, disk_region)

    # Calculate reaction rate.
    # TODO Include reaction probability?
    R = sigma * dv  # NOTE Simplification was made here.
    return R
