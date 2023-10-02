from pathlib import Path
from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, DiscreteRadialAxis
from collision import collision_outcome_probabilities
from collision import collision_rate
from config import Config
from disk import Disk, DiskRegion
from dust.relative_velocity import relative_velocity
from utils.functions import heaviside_theta


class Kernel():
    __slots__ = [
        "cfg", "mg", 
                  "R_coag",      "R_frag",
        "K",      "K_coag",      "K_frag", 
        "K_gain", "K_coag_gain", "K_frag_gain", 
        "K_loss", "K_coag_loss", "K_frag_loss",
    ]

    def __init__(
        self, 
        cfg: Config, 
        ijs: Optional[list[tuple[int, int]]] = None,
        R_coag: Optional[np.ndarray] = None,
        R_frag: Optional[np.ndarray] = None,
    ):
        self.cfg = cfg

        # Define discrete axes...
        rg = DiscreteRadialAxis(cfg) # ...for distance from central star,
        mg = DiscreteMassAxis(cfg)   # ...for particle mass,
        mc = mg.bin_centers
        self.mg = mg

        # TODO Correct? (only needed for log?)
        if cfg.enable_cancellation_handling and cfg.mass_axis_scale == "log":
            assert (mc[1:] / mc[:-1]).max() < 2, "Error: mass grid too coarse."

        # If relevant particle pairs are not specified explicitly,
        # assume all of them have to be taken into account.
        if ijs is None:
            ijs = []
            for i in range(mg.N):
                for j in range(mg.N):
                    ijs.append((i, j))

        if R_coag is R_frag is None:
            # Define PPD, & radial position of interest in it.
            disk = Disk(cfg, rg, mg) 
            disk_region = DiskRegion(cfg, disk)
            # Define relative dust particle velocities, & collision rates.
            dv = relative_velocity(cfg, disk, disk_region)
            R_coll = collision_rate(cfg, disk, disk_region)
            # Define probabilities for coagulation & fragmentation events.
            P_coag, P_frag = collision_outcome_probabilities(cfg, dv)
            # Define rate of coag./frag. events.
            R_coag, R_frag = R_coll * P_coag, R_coll * P_frag
        self.R_coag, self.R_frag = R_coag, R_frag

        # Initialize kernel matrices with zeros.
        self.K_coag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_coag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_frag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_frag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_gain = np.zeros(shape=[mg.N] * 3)
        self.K_loss = np.zeros(shape=[mg.N] * 3)
        self.K = np.zeros(shape=[mg.N] * 3)

        # Define gain & loss kernel sub-components for...
        # ...stick-and-hit coagulation processes.
        if cfg.enable_coagulation:
            K_coag_gain, K_coag_loss = self._K_coag(ijs)  # TODO Rename: K -> X?
            K_coag_gain *= R_coag
            K_coag_loss *= R_coag
            self.K_coag_gain += K_coag_gain
            self.K_coag_loss += K_coag_loss
            self.K_gain += K_coag_gain
            self.K_loss += K_coag_loss
        # ...fragmentation processes.
        if cfg.enable_fragmentation:
            K_frag_gain, K_frag_loss = self._K_frag(ijs)
            K_frag_gain *= R_frag
            K_frag_loss *= R_frag
            self.K_frag_gain += K_frag_gain
            self.K_frag_loss += K_frag_loss
            self.K_gain += K_frag_gain
            self.K_loss += K_frag_loss
        # Define total coagulation & fragmentation sub-kernels.
        self.K_coag = P_coag * self.K_coag_gain + P_coag * self.K_coag_loss
        self.K_frag = P_frag * self.K_frag_gain + P_frag * self.K_frag_loss
        # Define total kernel
        self.K += P_coag * self.K_coag
        self.K += P_frag * self.K_frag

    def _K_coag(
        self, 
        ijs: list[tuple[int, int]],
    ):
        mg = self.mg
        mc, N_m = mg.bin_centers, mg.N

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.
        for i, j in ijs:
            m_i, m_j = mc[i], mc[j]

            th = heaviside_theta(i - j)

            # Calculate combined mass after hit-and-stick collision.
            m_k = m_i + m_j

            # If a non-linear grid is used, the corresponding index will
            # not necessarily be an integer. Therefore, the resulting mass
            # has to be split onto the two neighboring bins, with indices:
            k_l = mg.index_from_value(m_k)  # := index of next-lower bin
            k_h = k_l + 1                   # := index of next-higher bin
            if k_h >= N_m:
                continue

            # Calculate masses corresponding to the two neighboring bins.
            m_l = mg.value_from_index(k_l)
            m_h = mg.value_from_index(k_h)
            assert m_k >= m_l and m_k <= m_h

            # Decide whether near-zero cancellation handling is required.
            might_cancel = (k_l == i)
            handle_cancellation = self.cfg.enable_cancellation_handling and might_cancel

            # Calculate fraction of mass "overflowing" into bin `k_h`.
            if handle_cancellation:
                eps = m_j / (m_h - m_l)  # Subtract analytically.
            else:
                eps = (m_i + m_j - m_l) / (m_h - m_l)

            # Check whether one of the masses is in the upper-most bin.
            near_upper_bound = i >= N_m - 1 or j >= N_m - 1

            # Subtract "loss" term from kernel.
            # ─────────────────────────────────────────────────────────────
            # Handle upper mass grid boundary.
            if not near_upper_bound:
                # Handle cancellation.
                if handle_cancellation:
                    K_loss[k_l, i, j] -= th * eps
                    K_loss[k_l, i, j] -= th if i == j else 0
                    # ^ TODO Why is this term here?
                    #        If removed, the solver crashes.
                # Handle "trivial" (non-cancelling) case.
                else:
                    K_loss[i, i, j] -= 1 if i < N_m - 1 else 0

            # Add "gain" term to kernel.
            # ─────────────────────────────────────────────────────────────
            # Handle upper mass grid boundary.
            if not near_upper_bound:
                # Handle cancellation.
                if handle_cancellation:
                    K_gain[k_h, i, j] += th * eps
                # Handle "trivial" (non-cancelling) case.
                else:
                    K_gain[k_l, i, j] += th * (1 - eps)
                    K_gain[k_h, i, j] += th * eps

        return K_gain, K_loss

    def _K_frag(
        self, 
        ijs: list[tuple[int, int]],
    ):
        mg = self.mg
        mc, N_m = mg.bin_centers, mg.N

        fragmentation_variant = self.cfg.fragmentation_variant

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        for i, j in ijs:
            m_i, m_j = mc[i], mc[j]
            th = heaviside_theta(i - j)

            # 1. Most basic/naive fragmentation implementation: "Pulverization"
            # ═════════════════════════════════════════════════════════════════

            if fragmentation_variant == "naive/pulverization":

                X = 0  # TODO Play around with this value, observe changes!
                k = 0
                m_k = mc[k]
                if min(i, j) > X:
                    eps = (m_i + m_j) / m_k
                    K_loss[i, i, j] -= 1
                    K_gain[k, i, j] += th * eps

            # 2. Mass redistribution following the MRN model
            # ═════════════════════════════════════════════════════════════════

            elif fragmentation_variant == "mrn":
                q = -11 / 6

                # Define total mass that needs to be "moved".
                m_tot = m_i + m_j

                # Define mass range resulting from fragmentation event.
                k_min = 0
                k_max = mg.index_from_value(m_tot)
                # k_max = mg.index_from_value(max(m_i, m_j))
                # ^ NOTE: This is a somewhat arbitrary choice:
                #   - One could also choose e.g. `m_max = max(m_i, m_j)`,
                #   - or something completely different, as long as mass is conserved.
                # ^ NOTE:
                #    - If we set `k_max = 1`, we expect the same behavior as 
                #      in "naive/pulverization" (with X set to ~= 1 there).

                # Calculate corresponding mass values from bin indices.
                m_min, m_max = mc[k_min], mc[k_max]
                assert m_min > mg.x_min
                assert m_max < mg.x_max

                # Calculate normalization constant for MRN distribution.
                S = sum([mc[k]**q for k in range(k_min, k_max)])  # TODO -> `k_max + 1` ? (below as well)
                assert S != 0  # TODO really needed?

                # Add mass to bins "receiving" mass in fragmentation event.
                for k in range(k_min, k_max):
                    A = mc[k]**q / S
                    K_gain[k, i, j] += m_tot / mc[k] * A * th
                    # K_gain[k, i, j] += m_tot * mc[k]**(q-1) / S * th
                    # ^ TODO: Why does this lead to changes in the kernel mass conservation plot?

                # Remove mass from bins corresponding to initial masses.
                K_loss[i, i, j] -= 1 

            else: 
                pass

        return K_gain, K_loss

    def save_to_file(
        self, 
        path: Optional[Path]=None
    ):
        K = self.K
        N = self.mg.N
        K = K.reshape((N, N*N))
        assert self.K.all() == K.reshape((N, N, N)).all()
        np.savetxt(str(path), K)
