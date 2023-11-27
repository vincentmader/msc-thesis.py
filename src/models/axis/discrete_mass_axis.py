from functions.physics.dust.dust_particle import particle_radius_from_mass

from .discrete_axis import DiscreteAxis


class DiscreteMassAxis(DiscreteAxis):
    __slots__ = ["particle_radii"]

    def __init__(self, cfg):
        super().__init__(
            cfg.mass_min_value,
            cfg.mass_max_value,
            cfg.mass_resolution,
            cfg.mass_axis_scale,
        )

        # Define field `particle_radii`.
        mc = self.bin_centers
        rho_s = cfg.dust_particle_density
        self.particle_radii = particle_radius_from_mass(mc, rho_s)
