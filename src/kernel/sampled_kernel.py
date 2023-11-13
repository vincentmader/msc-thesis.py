from typing import Optional

import numpy as np

from axis import DiscreteMassAxis
from config import Config
from kernel import Kernel


ALMOST_BUT_NOT_QUITE_ZERO = 1e-100  # TODO Choose this value. 


class SampledKernel(Kernel):
    __slots__ = [
        "P_ij",  #  = Probability of randomly selecting a collision pair $(i,j)$.
        "N_ij",  #  = Nr. of times that the collision $(i,j)$ was selected.
    ]

    def __init__(
        self, 
        cfg: Config, 
        N: np.ndarray, 
        W_ij: Optional[np.ndarray] = None,
        *args, **kwargs
    ):
        mg = DiscreteMassAxis(cfg)
        mc = mg.bin_centers
        if W_ij is None:  # Define weights, if not received as argument.
            kernel = Kernel(cfg)
            K = kernel.K
            W_ij = np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)])
            # TODO Use quadratic addition instead? (+ sqrt afterwards)

        # assert (np.abs(N[N < 0]) <= 1e-16).all(), N  # TODO Uncomment this line.
        N_i = np.abs(N[:, None])  # TODO `N` can be less than 0. Can I fix this?
        N_j = np.abs(N[None, :])

        m_i = mc[:, None]
        m_j = mc[None, :]

        P_ij = W_ij * N_i * N_j * m_i * m_j
        P_ij[P_ij == 0] = ALMOST_BUT_NOT_QUITE_ZERO
        P_ij = P_ij / P_ij.sum()  # Normalize. 
        assert np.abs(P_ij.sum() - 1) <= 1e-6

        self.P_ij = P_ij
        self.N_ij = np.zeros(shape=[N.shape[0]]*2)

        ijs = self._sample_ijs(cfg)
        super().__init__(cfg, ijs=ijs, *args, **kwargs)

    def _sample_ijs(self, cfg: Config) -> list[tuple[int, int]]:
        P_ij = self.P_ij
    
        assert len(P_ij.shape) == 2
        N_i, N_j = P_ij.shape
        assert N_i == N_j
        indices = range(N_i * N_j)
        P_ij = P_ij.reshape(N_i * N_j)
    
        assert (P_ij != 0).all()

        # If not sampling over all collisions, exclude "irrelevant" collisions $(i,j)$.
        # The "relevant" collisions are those with a probability significantly higher than 1e-100.
        N_relevant = np.sum(P_ij > ALMOST_BUT_NOT_QUITE_ZERO)
        N_relevant = np.sum(P_ij > 1e-20)
        N_sample = min(cfg.nr_of_samples, N_relevant)
        # print(N_relevant, N_sample)
        # N_sample = cfg.nr_of_samples
        # TODO Correctly do this. Uncommenting the lines above leads to `K_sampled != K_unsampled`.

        sampled = np.random.choice(indices, p=P_ij, size=N_sample, replace=False)
    
        ijs = []
        for s in sampled:
            i, j = s // N_i, s % N_i
            ijs.append((i, j))
            self.N_ij[i, j] += 1
        return ijs
