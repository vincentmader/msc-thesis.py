from typing import Optional

import numpy as np

from axis import DiscreteMassAxis, DiscreteRadialAxis
from collision import collision_outcome_probabilities_from_cutoff_velocity
from collision import collision_outcome_probabilities_from_maxwell_boltzmann
from collision import collision_rate
from config import Config
from disk import Disk, DiskRegion
from dust.relative_velocity import relative_velocity
from utils.functions import heaviside_theta


class Kernel():
    __slots__ = [
        "cfg", "mg", "K", "K_gain", "K_loss",
        "K_coag", "K_coag_gain", "K_coag_loss",
        "K_frag", "K_frag_gain", "K_frag_loss",
    ]

    def __init__(
        self, 
        cfg: Config, 
        ijs: Optional[list[tuple[int, int]]] = None,
    ):
        self.cfg = cfg

        # Define discrete axes for...
        # ...radial distance of disk region of interest from central star.
        rg = DiscreteRadialAxis(cfg)
        # ...particle mass.
        mg = DiscreteMassAxis(cfg)
        mc = mg.bin_centers
        self.mg = mg

        # If relevant particle pairs are not specified explicitly,
        # assume all of them have to be taken into account.
        if ijs is None:
            ijs = []
            for i in range(mg.N):
                for j in range(mg.N):
                    ijs.append((i, j))

        # Define protoplanetary disk, & the position of interest in it.
        disk = Disk(cfg, rg, mg)
        disk_region = DiskRegion(cfg, disk)

        # Define total relative dust particle velocities.
        dv = relative_velocity(cfg, disk, disk_region)

        # Define particle collision rates.
        if cfg.enable_physical_collisions:
            R_coll = collision_rate(cfg, disk, disk_region)
        else:  # In the most simple case, the rates are just set to 1.
            R_coll = np.ones(shape=[mg.N] * 2)

        if cfg.enable_cancellation_handling and cfg.mass_axis_scale == "log":
            # TODO Correct? (only needed for log?)
            assert (mc[1:] / mc[:-1]).max() < 2, "Error: mass grid too coarse."

        # Define probabilities for coagulation & fragmentation events.
        # Case 1: Include only hit-and-stick coagulation.
        if (cfg.enable_coagulation) and (not cfg.enable_fragmentation):
            P_coag, P_frag = 1, 0
        # Case 2: Include only fragmentation.
        elif (not cfg.enable_coagulation) and (cfg.enable_fragmentation):
            P_coag, P_frag = 0, 1
        # Case 3: Include neither (useless, but included for completeness).
        elif (not cfg.enable_coagulation) and (not cfg.enable_fragmentation):
            P_coag, P_frag = 0, 0
        # Case 4: Include both.
        else:
            # Case 4.1: Assume fragmentation if `v > v_cutoff`, else coagulation.
            if cfg.collision_outcome_variant == "cutoff_velocity":
                P_coag, P_frag = collision_outcome_probabilities_from_cutoff_velocity(cfg, dv)
            # Case 4.2: Calculate outcome probabilities from cutoff velocity & M.B. distribution.
            elif cfg.collision_outcome_variant == "mb_dist":
                P_coag, P_frag = collision_outcome_probabilities_from_maxwell_boltzmann(cfg, dv)
            # Case 4.3: Assume equal probabilities.
            else:
                P_coag, P_frag = 0.5, 0.5

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
            K_coag = self._K_coag(R_coll, ijs)
            self.K_coag_gain += P_coag * K_coag["gain"]
            self.K_coag_loss += P_coag * K_coag["loss"]
            self.K_gain += P_coag * K_coag["gain"]
            self.K_loss += P_coag * K_coag["loss"]
        # ...fragmentation processes.
        if cfg.enable_fragmentation:
            K_frag = self._K_frag(R_coll, ijs)
            self.K_frag_gain += P_frag * K_frag["gain"]
            self.K_frag_loss += P_frag * K_frag["loss"]
            self.K_gain += P_frag * K_frag["gain"]
            self.K_loss += P_frag * K_frag["loss"]
        # Define total coagulation & fragmentation sub-kernels.
        self.K_coag = P_coag * self.K_coag_gain + P_coag * self.K_coag_loss
        self.K_frag = P_frag * self.K_frag_gain + P_frag * self.K_frag_loss
        # Define total kernel
        self.K += P_coag * self.K_coag
        self.K += P_frag * self.K_frag

    def _K_coag(
        self, 
        R_coll: np.ndarray, 
        ijs: list[tuple[int, int]],
    ):

        mg = self.mg
        mc = mg.bin_centers
        m_max = mg.x_max  # Note: This is NOT the same as `mc[-1]`.
        N_m = mg.N

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.
        for i, j in ijs:
            m_i, m_j = mc[i], mc[j]

            th = heaviside_theta(i - j)

            # Get collision rate for mass pair.
            R = R_coll[i, j]

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
                    K_loss[k_l, i, j] -= R * th * eps
                    K_loss[k_l, i, j] -= R * th if i == j else 0
                    # ^ TODO Why is this term here?
                    #        If removed, the solver crashes.
                # Handle "trivial" (non-cancelling) case.
                else:
                    K_loss[i, i, j] -= R if i < N_m - 1 else 0

            # Add "gain" term to kernel.
            # ─────────────────────────────────────────────────────────────
            # Handle upper mass grid boundary.
            if not near_upper_bound:
                # Handle cancellation.
                if handle_cancellation:
                    K_gain[k_h, i, j] += R * th * eps
                # Handle "trivial" (non-cancelling) case.
                else:
                    K_gain[k_l, i, j] += R * th * (1 - eps)
                    K_gain[k_h, i, j] += R * th * eps

        return {"gain": K_gain, "loss": K_loss}

    def _K_frag(
        self, 
        R_coll: np.ndarray,
        ijs: list[tuple[int, int]],
    ):

        mg = self.mg
        mc = mg.bin_centers
        N_m = mg.N

        fragmentation_variant = self.cfg.fragmentation_variant

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        for i, j in ijs:
            m_i, m_j = mc[i], mc[j]

            th = heaviside_theta(i - j)
            
            R = R_coll[i, j]

            # 1. Most basic/naive fragmentation implementation: "Pulverization"
            # ═════════════════════════════════════════════════════════════════

            if fragmentation_variant == "naive/pulverization":

                X = 0  # TODO Play around with this value, observe changes!
                k = 0
                m_k = mc[k]
                if min(i, j) > X:
                    eps = (m_i + m_j) / m_k
                    K_loss[i, i, j] -= R
                    K_gain[k, i, j] += R * th * eps

            # 2. Mass redistribution following the MRN model
            # ═════════════════════════════════════════════════════════════════

            elif fragmentation_variant == "mrn":
                q = -11 / 6

                # Define total mass that needs to be "moved".
                m_tot = m_i + m_j

                # Define mass range resulting from fragmentation event.
                k_min = 0
                k_max = mg.index_from_value(m_tot)
                # ^ NOTE: This is a somewhat arbitrary choice:
                #   - One could also choose e.g. `m_max = max(m_i, m_j)`,
                #   - or something completely different, as long as mass is conserved.
                # ^ NOTE:
                #    - If we set `k_max = 1`, we expect the same behavior as 
                #      in "naive/pulverization" (with X set to ~= 1 there).

                # Calculate corresponding mass values from bin indices.
                m_min = mc[k_min]
                m_max = mc[k_max]
                assert m_min > mg.x_min
                assert m_max < mg.x_max

                # Calculate normalization constant for MRN distribution.
                S = sum([mc[k]**q for k in range(k_min, k_max)])
                assert S != 0

                # Add mass to bins "receiving" mass in fragmentation event.
                for k in range(k_min, k_max):
                    m_k = mc[k]
                    A = m_k**q / S
                    K_gain[k, i, j] += R * m_tot / m_k * A * th

                # Remove mass from bins corresponding to initial masses.
                K_loss[i, i, j] -= R 

            else: 
                pass

        return {"gain": K_gain, "loss": K_loss}
