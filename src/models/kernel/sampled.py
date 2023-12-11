from dataclasses import dataclass
from typing import Optional

import numpy as np

from config import Config
from models.axis import DiscreteMassAxis
from models.kernel import Kernel


ALMOST_BUT_NOT_QUITE_ZERO = 1e-100  # TODO Choose this value. 


@dataclass
class SampledKernel(Kernel):

    P_ij: np.ndarray  # = probability of randomly selecting a collision pair $(i,j)$.
    N_ij: np.ndarray  # = nr. of times that the collision $(i,j)$ was selected.
    W_ij: np.ndarray
    N_sample: int     # = nr. of kernel entries to sample
                      #   NOTE: This is NOT to be mistaken for `N_ij`.

    def __init__(
        self, 
        cfg:        Config, 
        N:          np.ndarray, 
        R_coag:     Optional[np.ndarray] = None,
        R_frag:     Optional[np.ndarray] = None,
        W_ij:       Optional[np.ndarray] = None,
        N_sample:   Optional[int]        = None,
        *args, 
        **kwargs
    ):
        # Define mass axis.
        mg = DiscreteMassAxis(cfg)
        mc = mg.bin_centers

        # Define weights, if not received as argument.
        if W_ij is None:
            K = Kernel(cfg, R_coag=R_coag, R_frag=R_frag).K
            self.W_ij = np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)])
            # W_ij = np.sum([mc[k] * K[k]**2 for k in range(mg.N)])**.5  # TODO Use lin. or quad. addition?
        else:
            self.W_ij = W_ij

        # Define nr. of particles per bin (& unit volume).
        m_i, m_j =  mc[:, None], mc[None, :]
        N_i, N_j =  N[:, None],   N[None, :]
        assert N.all() >= 0

        # Define sampling probability distribution.
        P_ij = self.W_ij * N_i * N_j * m_i * m_j
        # P_ij = W_ij * N_i**2 * N_j**2 * m_i * m_j

        # Exclude the lower-right (here "upper"?) matrix half from sampling.
        P_ij[np.triu_indices(mg.N, k=+1)] = 0

        # If sampling over all collisions, make sure that `P_ij  > 0` for all `i,j`.
        P_ij[P_ij == 0] = ALMOST_BUT_NOT_QUITE_ZERO

        # Normalize probability distribution.
        P_ij = P_ij / P_ij.sum()  
        assert np.abs(P_ij.sum() - 1) <= 1e-6

        self.P_ij = P_ij
        self.N_ij = np.zeros(shape=[N.shape[0]]*2)
        self.N_sample = cfg.nr_of_samples if N_sample is None else N_sample

        ijs = self._sample_ijs(cfg)
        super().__init__(cfg, R_coag, R_frag, ijs=ijs, *args, **kwargs)

    def _sample_ijs(self, cfg: Config) -> list[tuple[int, int]]:
        P_ij = self.P_ij
    
        assert len(P_ij.shape) == 2
        N_i, N_j = P_ij.shape
        assert N_i == N_j
        indices = range(N_i * N_j)
        P_ij = P_ij.reshape(N_i * N_j)

        # If sampling over all collisions, make sure that probability is > 0 everywhere.
        if self.N_sample == cfg.mass_resolution**2:
            assert (P_ij != 0).all()
            self.N_sample = cfg.nr_of_samples
        # If not sampling over all collisions, exclude "irrelevant" collisions $(i,j)$.
        # The "relevant" collisions are those with a probability significantly higher than 1e-100.
        else:
            N_relevant = np.sum(P_ij > 1e-16)
            self.N_sample = min(self.N_sample, N_relevant)

        assert P_ij.all() > 0
        sampled = np.random.choice(indices, p=P_ij, size=self.N_sample, replace=False)
    
        ijs = []
        for s in sampled:
            i = s // N_i
            j = s %  N_i
            ijs.append((i, j))
            self.N_ij[i, j] += 1
        return ijs
