import numpy as np

from disk import MassGrid, Disk, RadialGrid, DiskRegion
from dust.collision_rate import collision_rate
from utils.functions import heaviside_theta


class Kernel():

    def __init__(self, cfg):

        self.cfg = cfg

        # Define discrete axes for...
        # ...radial distance of disk region of interest from central star.
        rg = RadialGrid(cfg)
        # ...particle mass.
        mg = MassGrid(cfg)
        # Define protoplanetary disk, & the position of interest in it.
        disk = Disk(cfg, rg, mg)
        disk_region = DiskRegion(cfg, disk)
        # Define particle collision rates.
        if cfg.enable_physical_relative_velocities:
            R_coll = collision_rate(cfg, disk, disk_region)
        else:  # In the most simple case, the rates are just set to 1.
            R_coll = np.ones(shape=[mg.N]*2)

        # Define gain & loss kernel sub-components for...
        # ...stick-and-hit coagulation processes.
        K_coag = self._K_coag(mg, R_coll)
        K_coag_gain = K_coag["gain"]
        K_coag_loss = K_coag["loss"]
        K_coag = K_coag_gain + K_coag_loss
        # ...fragmentation processes.
        K_frag = self._K_frag(mg, R_coll)
        K_frag_gain = K_frag["gain"]
        K_frag_loss = K_frag["loss"]
        K_frag = K_frag_gain + K_frag_loss
        # Define total kernel.
        K = np.zeros(shape=[mg.N] * 3)
        K += K_coag if cfg.enable_coagulation else 0
        K += K_frag if cfg.enable_fragmentation else 0

        # Save kernel components to class object fields.
        self.K_coag_gain = K_coag_gain
        self.K_coag_loss = K_coag_loss
        self.K_frag_gain = K_frag_gain
        self.K_frag_loss = K_frag_loss
        self.K_coag = K_coag
        self.K_frag = K_frag
        self.K = K

    def _K_coag(self, mg, R_coll):

        N_m = mg.N
        m_max = mg.value_from_index(N_m - 1)

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.
        masses = mg.grid_cell_boundaries()[:-1]
        for i, m_i in enumerate(masses):
            for j, m_j in enumerate(masses):
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
                k_l = mg.index_from_value(m_k)  #  := index of next-lower bin
                k_h = k_l + 1                   #  := index of next-higher bin

                # Calculate masses corresponding to the two neighboring bins.
                m_l = mg.value_from_index(k_l)
                m_h = mg.value_from_index(k_h)

                # Decide whether near-zero cancellation handling is required.
                might_cancel = (k_l == i)
                handle_cancellation = self.cfg.enable_cancellation_handling and might_cancel

                # Calculate fraction of mass "overflowing" into bin `k_h`.
                if handle_cancellation:
                    eps = m_j / (m_h - m_l)  # Subtract analytically.
                else:
                    eps = (m_i + m_j - m_l) / (m_h - m_l)

                # Check whether one of the masses is in the upper-most 2 bins.
                near_upper_bound = i>=N_m-2 or j>=N_m-2

                # Subtract "loss" term from kernel.
                # ─────────────────────────────────────────────────────────────
                # Handle upper mass grid boundary.
                if not near_upper_bound:
                    # Handle cancellation.
                    if handle_cancellation:
                        K_loss[k_l, i, j] -= R * th * eps
                        K_loss[k_l, i, j] -= R/2 if i==j else 0
                    # Handle "trivial" (non-cancelling) case.
                    else:
                        K_loss[i, i, j] -= R if i<N_m-2 else 0

                # Add "gain" term to kernel.
                # ─────────────────────────────────────────────────────────────
                # Handle upper mass grid boundary.
                k_l = N_m-2 if k_l > N_m-2 else k_l
                k_h = N_m-2 if k_h > N_m-2 else k_h
                if not near_upper_bound:
                    # Handle cancellation.
                    if handle_cancellation:
                        K_gain[k_h, i, j] += R * th * eps
                    # Handle "trivial" (non-cancelling) case.
                    else:
                        K_gain[k_l, i, j] += R * th * (1 - eps)
                        K_gain[k_h, i, j] += R * th * eps

        return {"gain": K_gain, "loss": K_loss}

    def _K_frag(self, mg, R_coll):

        masses = mg.grid_cell_centers()
        boundaries = mg.grid_cell_boundaries()
        dm = boundaries[1:] - boundaries[:-1]
        N_m = mg.N

        fragmentation_variants = self.cfg.enable_fragmentation_variant

        K_gain = np.zeros(shape=[N_m] * 3)
        K_loss = np.zeros(shape=[N_m] * 3)

        for i, m_i in enumerate(masses):
            for j, m_j in enumerate(masses):
                th = heaviside_theta(i - j)

                # 1. Most basic/naive fragmentation implementation: "Pulverization"
                # ═════════════════════════════════════════════════════════════════

                if "naive/pulverization" in fragmentation_variants:

                    X = 10 # TODO Play around with this value, observe changes!
                    k = 0
                    m_k = masses[k]
                    if min(i, j) > X:
                        eps = (m_i + m_j) / m_k
                        K_loss[i, i, j] -= R_coll[i, j]
                        K_gain[k, i, j] += R_coll[i, j] * th * eps

                # 2. Mass redistribution following the MRN model
                # ═════════════════════════════════════════════════════════════════

                if "mrn" in fragmentation_variants:

                    q = -11/6

                    # Define total mass that needs to be "moved".
                    m_tot = m_i + m_j

                    # Define range of masses resulting from fragmentation event.
                    m_min = masses[0]
                    m_max = m_tot
                    k_min = mg.index_from_value(m_min)
                    k_max = mg.index_from_value(m_max)
                    k_max = min(N_m-1, k_max)

                    # This is the lowest bin at which fragmentation can occur.
                    if min(i, j) < k_min:
                        continue

                    top = m_i * dm[i] * R_coll[i, j] + m_j * dm[j] * R_coll[i, j]
                    bottom = sum([masses[k] * dm[k] * masses[k]**q for k in range(k_min, k_max)])
                    A = top / bottom

                    # Remove mass from bins corresponding to initial masses.
                    K_loss[i, i, j] -= R_coll[i, j]

                    # Add mass to bins corresponding to resulting masses.
                    # Loop over all bins that are "receiving" mass.
                    for k, m_k in enumerate(masses):
                        if m_k > m_max:
                            continue

                        # Add mass to bin.
                        n = A * m_k**q
                        K_gain[k, i, j] += R_coll[i, j] * n

        return {"gain": K_gain, "loss": K_loss}
