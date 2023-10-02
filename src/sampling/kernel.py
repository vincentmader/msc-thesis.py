from typing import Optional

import numpy as np

from config import Config
from kernel import Kernel


class SampledKernel(Kernel):
    __slots__ = ["P_ij"]

    def __init__(
        self, 
        cfg: Config, 
        N: np.ndarray, 
        W_ij: Optional[np.ndarray]=None,
        *args, **kwargs
    ):

        if W_ij is None:
            kernel = Kernel(cfg)
            K, mg = kernel.K, kernel.mg
            mc = mg.bin_centers
            W_ij = sum([mc[k] * np.abs(K[k]) for k in range(mg.N)])

        N_i = np.abs(N[:, None])
        N_j = np.abs(N[None, :])

        P_ij = W_ij * N_i * N_j  # TODO Is this multiplication correct?
        P_ij = P_ij / P_ij.sum()
        assert np.abs(P_ij.sum() - 1) <= 1e-6  # TODO -> 1e-16 ?
        self.P_ij = P_ij
        
        ijs = self._sample_ijs(cfg)

        super().__init__(cfg, ijs=ijs, *args, **kwargs)

    def _sample_ijs(self, cfg: Config) -> list[tuple[int, int]]:
        P_ij = self.P_ij
    
        assert len(P_ij.shape) == 2
        N_i, N_j = P_ij.shape
        assert N_i == N_j
        indices = range(N_i * N_j)
        P_ij = P_ij.reshape(N_i * N_j)
    
        N_sample = cfg.nr_of_samples
        sampled = np.random.choice(
            indices, p=P_ij, size=N_sample, # NOTE: Pairs can be selected multiple times
        )
    
        ijs = []
        for s in sampled:
            i = s // N_i
            j = s % N_i
            if (i, j) in ijs: # NOTE: This helped a lot! (no duplicates)
                continue
            ijs.append((i, j))
        return ijs
