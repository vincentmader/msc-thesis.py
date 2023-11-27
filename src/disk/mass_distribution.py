import numpy as np

from models.axis import DiscreteMassAxis, DiscreteRadialAxis
from disk import Disk


def dirac_delta(cfg):
    mg = DiscreteMassAxis(cfg)
    rg = DiscreteRadialAxis(cfg)
    disk = Disk(cfg, rg, mg)
    r = rg.bin_centers  
    rho_g = disk.midplane_gas_volume_density

    x_0 = cfg.initial_mass_bin

    m = mg.bin_centers
    dm = mg.bin_widths[x_0]

    R = cfg.distance_to_star
    if cfg.enable_physical_gas_density:
        RHO_g = np.interp(R, r, rho_g)
        # Is linear interpolation good enough here? -> Yes!
    else:
        RHO_g = 1

    M = RHO_g / m[x_0] / dm

    # COMPARISON: Kees' code
    # Create initial distribution: Nr of particles per interval dmass per volume.
    # n_dust = np.zeros(N_m)
    # n_dust[0] = 1.0 / dmgrain[0]

    n0 = np.zeros([mg.N])
    n0[x_0] = M

    return n0


def mrn_distribution(cfg):
    mg = DiscreteMassAxis(cfg)
    rg = DiscreteRadialAxis(cfg)
    disk = Disk(cfg, rg, mg)
    r = rg.bin_centers 
    rho_g = disk.midplane_gas_volume_density

    m = mg.bin_centers
    dm = mg.bin_widths[x_0]

    m_min = m[0]
    m_max = m[-1]

    R = cfg.distance_to_star
    RHO_g = np.interp(R, r, rho_g)

    M = RHO_g / m[x_0] / dm

    q = -11 / 6
    n0 = np.ones(mg.N) * m**q

    A = 1 / (q + 2) * (m_max**(q + 2) - m_min**(q + 2))
    n0 = n0 * M / A

    return n0
