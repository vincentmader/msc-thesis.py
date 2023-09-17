import numpy as np

from axis import DiscreteMassAxis
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel.kernel import KernelSubplot


class KernelMassConservationSubplot(KernelSubplot):
    __slots__ = ["sum_ij"]

    def __init__(
        self,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        *args, **kwargs
    ):
        self.sum_ij = test_mass_conservation(mg, K)
        super().__init__(mg, self.sum_ij, *args, **kwargs)
