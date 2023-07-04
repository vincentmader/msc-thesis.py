import numpy as np

from axis import DiscreteMassAxis, DiscreteRadialAxis
from disk import Disk, DiskRegion
from dust.collision_rate import collision_rate
from utils.functions import heaviside_theta


class Kernel():

    def __init__(self, cfg):

        self.cfg = cfg

        # Define discrete axes for...
        # ...radial distance of disk region of interest from central star.
        rg = DiscreteRadialAxis(cfg)
        # ...particle mass.
        mg = DiscreteMassAxis(cfg)
        mc = mg.grid_cell_centers
        self.mg = mg
        # Define protoplanetary disk, & the position of interest in it.
        disk = Disk(cfg, rg, mg)
        disk_region = DiskRegion(cfg, disk)
        # Define particle collision rates.
        if cfg.enable_physical_relative_velocities:
            R_coll = collision_rate(cfg, disk, disk_region)
        else:  # In the most simple case, the rates are just set to 1.
            R_coll = np.ones(shape=[mg.N] * 2)

        if cfg.enable_cancellation_handling and cfg.mass_axis_scale == "log":
            # TODO Correct? (only needed for log?)
            assert (mc[1:] / mc[:-1]).max() < 2, "Error: mass grid too coarse."

        # Initialize kernel matrices with zeros.
        self.K_coag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_coag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_frag_gain = np.zeros(shape=[mg.N] * 3)
        self.K_frag_loss = np.zeros(shape=[mg.N] * 3)
        self.K_gain = np.zeros(shape=[mg.N] * 3)
        self.K_loss = np.zeros(shape=[mg.N] * 3)

        # Define gain & loss kernel sub-components for...
        # ...stick-and-hit coagulation processes.
        K_coag = self._K_coag(R_coll)
        if cfg.enable_coagulation:
            self.K_coag_gain += K_coag["gain"]
            self.K_coag_loss += K_coag["loss"]
            self.K_gain += K_coag["gain"]
            self.K_loss += K_coag["loss"]
        # ...fragmentation processes.
        K_frag = self._K_frag(R_coll)
        if cfg.enable_fragmentation:
            self.K_frag_gain += K_frag["gain"]
            self.K_frag_loss += K_frag["loss"]
            self.K_gain += K_frag["gain"]
            self.K_loss += K_frag["loss"]
        # Define total kernel.
        self.K_coag = self.K_coag_gain + self.K_coag_loss
        self.K_frag = self.K_frag_gain + self.K_frag_loss
        self.K = self.K_coag + self.K_frag

    def _K_coag(self, R_coll):

        mg = self.mg
        mc = mg.grid_cell_centers
        m_max = mg.x_max  # Note: This is NOT the same as `mc[-1]`.
        N_m = mg.N

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.
        for i, m_i in enumerate(mc):
            for j, m_j in enumerate(mc):
                th = heaviside_theta(i - j)

                # Get collision rate for mass pair.
                R = R_coll[i, j]

                # Calculate combined mass after hit-and-stick collision.
                m_k = m_i + m_j
                if m_k >= m_max:
                    continue

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

                # Check whether one of the masses is in the upper-most 2 bins.
                near_upper_bound = i >= N_m - 2 or j >= N_m - 2

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
                        K_loss[i, i, j] -= R if i < N_m - 2 else 0

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

    def _K_frag(self, R_coll):

        mg = self.mg
        mc = mg.grid_cell_centers
        dm = mg.grid_cell_widths
        N_m = mg.N

        fragmentation_variants = self.cfg.enable_fragmentation_variant

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        for i, m_i in enumerate(mc):
            for j, m_j in enumerate(mc):
                th = heaviside_theta(i - j)
                
                R = R_coll[i, j]

                # 1. Most basic/naive fragmentation implementation: "Pulverization"
                # ═════════════════════════════════════════════════════════════════

                if "naive/pulverization" in fragmentation_variants:

                    X = 10  # TODO Play around with this value, observe changes!
                    k = 0
                    m_k = mc[k]
                    if min(i, j) > X:
                        eps = (m_i + m_j) / m_k
                        K_loss[i, i, j] -= R
                        K_gain[k, i, j] += R * th * eps

                # 2. Mass redistribution following the MRN model
                # ═════════════════════════════════════════════════════════════════
                q = -11 / 6

                if "mrn" in fragmentation_variants:

                    # Consider the very first case `i == j == 0`.
                    # - Here, the variable `bottom` will be set to zero
                    #   due to summation over an empty array.
                    # - This will lead to a division by zero, i.e. `A = infinity`.
                    # - Therefore, skip this case.
                    if i == 0 and j == 0:
                        continue
                        # ^ TODO:
                        #   - When else can this happen?
                        #   - Handle those cases!

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

                    # Make sure that all resulting masses lie inside the grid.
                    # k_max = min(N_m - 1, k_max)  # NOTE Interestingly, this is not needed.
                    assert k_max < N_m

                    # Calculate corresponding mass values from bin indices.
                    m_min = mc[k_min]
                    m_max = mc[k_max]
                    assert m_min > mg.x_min
                    assert m_max < mg.x_max

                    # VARIANT A: This leads to a quite big mass loss.
                    # ─────────────────────────────────────────────────────────────────────────────

                    # Remove mass from bins corresponding to initial masses.
                    # K_loss[i, i, j] -= R

                    # # TODO Calculate value to add to indices "receiving" mass.
                    # top = m_i * dm[i] * R + m_j * dm[j] * R
                    # bottom = sum(
                    #     [mc[k] * dm[k] * mc[k]**q for k in range(k_min, k_max)]
                    # )
                    # assert bottom != 0
                    # A = top / bottom

                    # # Add mass to bins corresponding to resulting masses.
                    # # Loop over all bins that are "receiving" mass.
                    # for k in range(k_min, k_max):
                    #     m_k = mc[k]

                    #     # Add mass to bin.
                    #     K_gain[k, i, j] += R * m_k**q * A

                    # NOTE: Here, several changes needed to be applied:
                    # - Remove `dm` terms from `top`.
                    # - Remove `R` terms from `top`, they are applied in `K_gain`.
                    # - Interestingly: Subtract `R` from `K_loss[j, i, j]` as well!
                    #   + I didn't have to to that in coagulation (?)
                    #   + Did I do that, but am not seeing where?
                    # - ...?

                    # VARIANT B: This was just for testing, it does the same as "naive" variant.
                    # ─────────────────────────────────────────────────────────────────────────────

                    # K_loss[i, i, j] -= R
                    # ^ NOTE: I need to subtract from both `i` and `j`!
                    #   - Why is that?
                    #   - I did not have to dod that in the "naive/pulverization" variant.

                    # k = 0
                    # m_k = mc[k]
                    # K_gain[k, i, j] += R * m_tot / m_k * th

                    # VARIANT C: Here, mass is conserved quite well! 
                    # ─────────────────────────────────────────────────────────────────────────────

                    # Remove mass from bins corresponding to initial masses.
                    K_loss[i, i, j] -= R 

                    S = 0
                    for k in range(k_min, k_max):
                        m_k = mc[k]
                        S += m_k**q

                    for k in range(k_min, k_max):
                        m_k = mc[k]
                        A = m_k**q / S
                        K_gain[k, i, j] += R * m_tot / m_k * A * th

        return {"gain": K_gain, "loss": K_loss}
