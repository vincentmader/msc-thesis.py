from typing import Optional

import numpy as np

from config import Config
from models.axis import DiscreteMassAxis
from models.kernel import Kernel

def sampling_weights(
    W_ij: Optional[np.ndarray],
    cfg, mg, R_coag, R_frag
) -> np.ndarray:
    if W_ij is None:
        mc = mg.bin_centers
        K  = Kernel(cfg, R_coag=R_coag, R_frag=R_frag).K
        return np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)], axis=0)
    else:
        return W_ij

ALMOST_BUT_NOT_QUITE_ZERO = 1e-100  # TODO Choose this value. 
def sampling_probability(W_ij, N, mg):
    mc = mg.bin_centers
    m_i, m_j =  mc[:, None], mc[None, :]
    N_i, N_j =   N[:, None],  N[None, :]

    P_ij = W_ij * N_i * N_j * m_i * m_j

    P_ij[np.triu_indices(mg.N, k=+1)] = 0
    P_ij[P_ij == 0] = ALMOST_BUT_NOT_QUITE_ZERO

    P_ij = P_ij / P_ij.sum()  
    return P_ij

class SampledKernelV2(Kernel):

    N_sample: int     # = nr. of kernel entries to sample
    W_ij: np.ndarray
    P_ij: np.ndarray  # = probability of randomly selecting a collision pair $(i,j)$.
    N_ij: np.ndarray  # = nr. of times that the collision $(i,j)$ was selected.

    def __init__(self, 
        cfg:        Config, 
        N:          np.ndarray, 
        N_sample:   int, 
        R_coag:     Optional[np.ndarray] = None,
        R_frag:     Optional[np.ndarray] = None,
        W_ij:       Optional[np.ndarray] = None,
    ):
        mg = DiscreteMassAxis(cfg)

        self.W_ij     = sampling_weights(W_ij, cfg, mg, R_coag, R_frag)
        self.P_ij     = sampling_probability(self.W_ij, N, mg)
        self.N_ij     = np.zeros(shape=[N.shape[0]]*2)
        self.N_sample = N_sample

        ijs = self._sample_ijs(cfg)
        super().__init__(cfg, R_coag, R_frag, ijs=ijs)
        if cfg.allow_duplicate_sampling:
            self.K = np.array([K_k / self.P_ij for K_k in self.K])

    def _sample_ijs(self, cfg: Config) -> list[tuple[int, int]]:
        P_ij = self.P_ij
    
        N_i, N_j = P_ij.shape
        indices = range(N_i * N_j)
        P_ij = P_ij.reshape(N_i * N_j)

        print(self.N_sample)
        if self.N_sample == cfg.mass_resolution**2:
            assert (P_ij != 0).all()
        else:
            N_relevant = np.sum(P_ij[P_ij > 1e-16])
            # self.N_sample = min(self.N_sample, N_relevant)

        assert P_ij.all() > 0
        replace = cfg.allow_duplicate_sampling
        sampled = np.random.choice(indices, p=P_ij, size=self.N_sample, replace=replace)
    
        ijs = []
        for s in sampled:
            i = s // N_i
            j = s %  N_i
            ijs.append((i, j))
            self.N_ij[i, j] += 1
        return ijs
