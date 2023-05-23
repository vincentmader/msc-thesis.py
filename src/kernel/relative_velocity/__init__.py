import numpy as np

from config import Config
from disk import Disk, DiskRegion
from utils.functions import root_mean_squared
from .azimuthal import dv_azimuthal
from .brownian_motion import dv_brownian_motion
from .differential_settling import dv_differential_settling
from .radial_drift import dv_radial_drift
from .turbulence import dv_turbulence


def relative_velocity(
    cfg: Config,
    disk: Disk,
    disk_region: DiskRegion,
):
    if cfg.enable_physical_relative_velocities == []:
        # return 35 * np.ones(shape=[disk.mass_axis.N_x] * 2)
        return np.ones(shape=[disk.mass_axis.N_x] * 2)

    dvs = []
    for dv_name in cfg.enable_physical_relative_velocities:

        if dv_name == "brownian_motion":
            dv_br = dv_brownian_motion(cfg, disk, disk_region)
            dvs.append(dv_br)

        elif dv_name == "differential_settling":  # note: Ignore for now.
            dv_ds = dv_differential_settling(cfg, disk, disk_region)
            dvs.append(dv_ds)

        elif dv_name == "radial_drift":
            dv_rd = dv_radial_drift(cfg, disk, disk_region)
            dvs.append(dv_rd)

        elif dv_name == "turbulence":
            dv_tu = dv_turbulence(cfg, disk, disk_region)
            dvs.append(dv_tu)

        elif dv_name == "azimuthal":
            dv_az = dv_azimuthal(cfg, disk, disk_region)
            dvs.append(dv_az)

    return root_mean_squared(dvs)
