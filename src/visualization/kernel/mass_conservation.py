from dust import particle_radius_from_mass
from kernel import Kernel
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel.pcolor_matrix_subplot import PcolorMatrixSubplot


class KernelMassConservationSubplot(PcolorMatrixSubplot):
    __slots__ = ["kernel"]

    def __init__(self, kernel: Kernel, *args, **kwargs):

        self.kernel = kernel
        cfg, mg, K = kernel.cfg, kernel.mg, kernel.K
        mc = mg.grid_cell_centers
        rho_s = cfg.dust_particle_density
        ac = particle_radius_from_mass(mc, rho_s)

        x, y, z = ac, ac, test_mass_conservation(cfg, mg, K)
        super().__init__(x, y, z, *args, **kwargs)
