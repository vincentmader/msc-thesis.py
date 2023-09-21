import numpy as np

from axis import DiscreteMassAxis
from kernel.mass_conservation import test_mass_conservation
from visualization.kernel import KernelSubplot


class KernelMassConservationSubplot(KernelSubplot):
    __slots__ = ["err_matrix", "err_total"]

    def __init__(
        self,
        mg: DiscreteMassAxis,
        K: np.ndarray,
        *args, **kwargs
    ):
        self.err_matrix, self.err_total = test_mass_conservation(mg, K)

        kwargs["title"] = kwargs["title"]\
            if "title" in kwargs.keys()\
            else r"kernel mass error $\Delta_{ij}=\sum_k m_k\cdot K_{kij}$, $\sum_{ij}|\Delta_ij|$ = " + f"{self.err_total}"

        super().__init__(
            mg, self.err_matrix, 
            *args, **kwargs
        )
