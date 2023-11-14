from pathlib import Path
from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, DiscreteRadialAxis
from collision import collision_outcome_probabilities, collision_rate
from config import Config
from disk import Disk, DiskRegion
from dust.relative_velocity import relative_velocity
from utils.functions import heaviside_theta


class Kernel():
    __slots__ = [
        "cfg",    "mg", 
        "K",      "K_coag",      "K_frag", 
        "K_gain", "K_coag_gain", "K_frag_gain", 
        "K_loss", "K_coag_loss", "K_frag_loss",
                  "R_coag",      "R_frag",
    ]

    def __init__(self, 
        cfg:    Config, 
        ijs:    Optional[list[tuple[int, int]]] = None,
        R_coag: Optional[np.ndarray]            = None,
        R_frag: Optional[np.ndarray]            = None,
    ):
        self.cfg = cfg

        # Define discrete axes...
        rg = DiscreteRadialAxis(cfg) # ...for distance from central star,
        mg = DiscreteMassAxis(cfg)   # ...for particle mass,
        mc = mg.bin_centers
        self.mg = mg

        if cfg.enable_cancellation_handling and cfg.mass_axis_scale == "log":
            assert (mc[1:] / mc[:-1]).max() < 2, "Error: mass grid too coarse."

        # If relevant particle pairs are not specified explicitly:
        # -> Assume all of them have to be taken into account.
        if ijs is None:
            ijs = []
            for i in range(mg.N):
                for j in range(mg.N):
                    ijs.append((i, j))
        for (i, j) in ijs:
            if (j, i) not in ijs:
                ijs.append((j, i))

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
        if (R_coag is None and R_frag is not None) or (R_frag is None and R_coag is not None):
            raise Exception("Unhandled case: Either both `R_coag` & `R_frag` have to be given, or neither.")
        self.R_coag, self.R_frag = R_coag, R_frag

        # Initialize kernel matrices with zeros.
        self.K_coag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_coag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_frag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_frag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_gain      = np.zeros(shape=[mg.N] * 3)
        self.K_loss      = np.zeros(shape=[mg.N] * 3)
        self.K           = np.zeros(shape=[mg.N] * 3)

        # Define gain & loss kernel sub-components for...
        # ...stick-and-hit coagulation processes.
        if cfg.enable_coagulation:
            K_coag_gain, K_coag_loss = self._K_coag(ijs)  # TODO Rename: K -> X?
            K_coag_gain      *= R_coag
            K_coag_loss      *= R_coag
            self.K_coag_gain += K_coag_gain
            self.K_coag_loss += K_coag_loss
            self.K_gain      += K_coag_gain
            self.K_loss      += K_coag_loss
        # ...fragmentation processes.
        if cfg.enable_fragmentation:
            K_frag_gain, K_frag_loss = self._K_frag(ijs)
            K_frag_gain      *= R_frag
            K_frag_loss      *= R_frag
            self.K_frag_gain += K_frag_gain
            self.K_frag_loss += K_frag_loss
            self.K_gain      += K_frag_gain
            self.K_loss      += K_frag_loss
        # Define total coagulation & fragmentation sub-kernels.
        self.K_coag = self.K_coag_gain + self.K_coag_loss
        self.K_frag = self.K_frag_gain + self.K_frag_loss
        # Define total kernel
        self.K += self.K_coag
        self.K += self.K_frag

    def _K_coag(self, ijs: list[tuple[int, int]]):
        mg = self.mg
        mc, N_m = mg.bin_centers, mg.N
        dm = mg.bin_widths

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.
        for i, j in ijs:  # TODO Handle cases where i < j. (lower right of kernel = 0!)
            ii, jj = (i, j) if i >= j else (j, i)
            assert ii >= jj
            # NOTE: The variables `ii` or `jj` are used for the case where we
            #       have to flip indices, i.e. do `(i,j) -> (j,i)`. This "flip" 
            #       is needed to assure that the kernel lives entirely in the 
            #       top-left half of the matrix (due to symmetry of problem, 
            #       collision "ij = ji"). When doing the flip, we cannot do a 
            #       simple tuple deconstruction like `(i,j) = (j, i)`, since the 
            #       indices are needed elsewhere in their "unflipped" initial state.
            m_i, m_j = mc[ii], mc[jj]
            # NOTE: Ignoring a possible index flip (like mentioned above) when defining
            #       the masses `m_i` and `m_i` is possible when they're only used for
            #       defining their sum `m_tot` (commutative/symmetric).
            #       In general, this is not necessarily the case, but here it's alright.

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
            m_l, m_h = mg.value_from_index(k_l), mg.value_from_index(k_h)
            assert (m_k >= m_l) and (m_k <= m_h)

            # Decide whether near-zero cancellation handling is required.
            might_cancel = (k_l == i)
            handle_cancellation = (self.cfg.enable_cancellation_handling and might_cancel)

            # Calculate fraction of mass "overflowing" into bin `k_h`.
            if not handle_cancellation:
                eps = (m_i + m_j - m_l) / (m_h - m_l)
            else:
                eps = m_j / (m_h - m_l)  # Subtract analytically.

            # Check whether one of the masses is in the upper-most bin.
            near_upper_bound = (i >= N_m - 1) or (j >= N_m - 1)

            # Subtract "loss" term from kernel.
            # ─────────────────────────────────────────────────────────────
            if not near_upper_bound:
                if handle_cancellation:
                    K_loss[k_l, ii, jj] -= th * eps
                    # K_loss[k_l, ii, jj] -= th if i == j else 0
                    # ^ TODO Why is this term here? (If removed, the solver crashes/crashed)
                else:  # Handle "trivial" (non-cancelling) case.
                    K_loss[i, ii, jj] -= 1 if i < N_m - 1 else 0

            # Add "gain" term to kernel.
            # ─────────────────────────────────────────────────────────────
            if not near_upper_bound:
                if handle_cancellation:
                    K_gain[k_h, ii, jj] += th * eps
                else:  # Handle "trivial" (non-cancelling) case.
                    K_gain[k_l, ii, jj] += th * (1 - eps)
                    K_gain[k_h, ii, jj] += th * eps

        return K_gain, K_loss

    def _K_frag(self, ijs: list[tuple[int, int]]):
        mg = self.mg
        mc, N_m = mg.bin_centers, mg.N
        dm = mg.bin_widths

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        for i, j in ijs:  # TODO Handle cases where i < j. (lower right of kernel = 0!)
            ii,  jj  = (i, j) if i >= j else (j, i)
            m_i, m_j = mc[i], mc[j]
            th = heaviside_theta(i - j)

            # 1. Most basic/naive fragmentation implementation: "Pulverization"
            # ═════════════════════════════════════════════════════════════════

            if self.cfg.fragmentation_variant == "naive/pulverization":

                X = 0  # TODO Play around with this value, observe changes!
                k = 0
                m_k = mc[k]
                if min(i, j) > X:
                    eps = (m_i + m_j) / m_k
                    K_loss[i, ii, jj] -= 1
                    K_gain[k, ii, jj] += th * eps

            # 2. Mass redistribution following the MRN model
            # ═════════════════════════════════════════════════════════════════

            elif self.cfg.fragmentation_variant == "mrn":
                q = -11 / 6

                # Define total mass that needs to be "moved".
                m_tot = m_i + m_j

                # Define mass range resulting from fragmentation event.
                k_min, k_max = 0, mg.index_from_value(m_tot)
                # ^ NOTE: This is a somewhat arbitrary choice: One could also choose
                #   - e.g. `m_max = max(m_i, m_j)`, or
                #   - something completely different, as long as mass is conserved.
                # ^ NOTE:
                #    - If we set `k_max = 1`, we expect the same behavior as 
                #      in "naive/pulverization" (with X set to ~= 1 there).

                # Calculate corresponding mass values from bin indices.
                m_min, m_max = mc[k_min], mc[k_max]
                assert m_min > mg.x_min
                assert m_max < mg.x_max

                # Calculate normalization constant for MRN distribution.
                S = sum([ mc[kk]**(q+1) * dm[kk] for kk in range(k_min, k_max) ])
                # TODO ^ Multiply with `dm[kk]` inside sum?
                assert S != 0  # TODO Really needed? Can `S==0`? Can I skip if it does?
                
                # Add mass to bins corresponding to newly created particles.
                for k in range(k_min, k_max):
                    K_gain[k, ii, jj] += m_tot * mc[k]**q / S * th * dm[k]

                # Remove mass from bins corresponding to initial particles.
                K_loss[i, ii, jj] -= 1

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
