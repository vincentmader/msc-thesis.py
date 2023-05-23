import numpy as np

from disk import MassGrid, Disk, RadialGrid


x_0 = 0


def dirac_delta(cfg):
    mg = MassGrid(cfg)
    rg = RadialGrid(cfg)
    disk = Disk(cfg, rg, mg)
    r = rg.grid_cell_centers() # TODO: Use centers or bounds?
    Sigma_g = disk.gas_surface_density(r)
    rho_g = disk.midplane_gas_volume_density(r[:-1], Sigma_g)

    m = mg.grid_cell_centers()
    dm = mg.grid_cell_width()

    R = cfg.distance_to_star
    if cfg.enable_physical_gas_density:
        RHO_g = np.interp(R, r[:-1], rho_g) 
        # Is linear interpolation good enough here? -> Yes!
    else:
        RHO_g = 1
    DM = dm * m[x_0]

    M = RHO_g / m[x_0] / DM

    n0 = np.zeros([mg.N_x])
    n0[x_0] = M

    return n0


def mrn_distribution(cfg):
    mg = MassGrid(cfg)
    rg = RadialGrid(cfg)
    disk = Disk(cfg, rg, mg)
    r = rg.grid_cell_centers() # TODO: Use centers or bounds?
    Sigma_g = disk.gas_surface_density(r)
    rho_g = disk.midplane_gas_volume_density(r[:-1], Sigma_g)

    m = mg.grid_cell_centers()
    dm = mg.grid_cell_width()

    m_min = m[0]
    m_max = m[-1]

    R = cfg.distance_to_star
    RHO_g = np.interp(R, r[:-1], rho_g) 
    DM = dm * m[x_0]
    M = RHO_g / m[x_0] / DM

    q = -11/6
    n0 = np.ones(mg.N_x) * m**q

    A = 1/(q+2) * (m_max**(q+2) - m_min**(q+2))
    n0 = n0 * M / A

    return n0
