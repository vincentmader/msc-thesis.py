from typing import Optional

from axis import KernelAxis
from dust import particle_mass_from_radius
from kernel import Kernel
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel import KernelSubplot


class KernelMassConservationSubplot(KernelSubplot):
    __slots__ = ["kernel", "err_matrix", "err_total"]

    def __init__(
        self,
        kernel: Kernel,
        axis_label_variant: Optional[KernelAxis] = KernelAxis.Radius, # TODO Rename? -> `axis_variant`
        *args, **kwargs
    ):
        self.kernel = kernel
        self.err_matrix, self.err_total = test_mass_conservation(kernel)

        kwargs["title"] = kwargs["title"]\
            if "title" in kwargs.keys()\
            else r"kernel mass error $\Delta_{ij}=\sum_k m_k\cdot K_{kij}$, $\sum_{ij}|\Delta_ij|$ = " + f"{self.err_total}"

        super().__init__(self.kernel.mg, self.err_matrix, axis_label_variant=axis_label_variant, *args, **kwargs)

    def format_coord(self, x, y):
        if self.axis_variant is KernelAxis.Radius:
            rho_s = self.kernel.cfg.dust_particle_density
            m_i = particle_mass_from_radius(y, rho_s)
            m_j = particle_mass_from_radius(x, rho_s)
            i = self.mg.index_from_value(m_i)
            j = self.mg.index_from_value(m_j)
        elif self.axis_variant is KernelAxis.Radius:
            i = self.mg.index_from_value(y)
            j = self.mg.index_from_value(x)
        else:  # -> `KernelAxis.Bin`
            i, j = int(y), int(x)

        return f"{i = }, {j = }, sum_k K_kij = {self.err_matrix[i, j]:.2}"
