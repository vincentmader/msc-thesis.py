import numpy as np

from utils.functions import heaviside_theta


class Kernel():

    def __init__(self, cfg):
        self.cfg = cfg

    def K(self, mg, R_coll):
        cfg = self.cfg

        N_m = mg.N_x
        # Create 3D array of 0s with `N_m` entries in each dimension.
        K = np.zeros(shape=[N_m] * 3)
        # Apply coagulation & fragmentation processes.
        K += self.K_coag(mg, R_coll) if cfg.enable_coagulation else 0
        K += self.K_frag(mg, R_coll) if cfg.enable_fragmentation else 0
        return K

    def K_coag(self, mg, R_coll):
        N_m = mg.N_x

        K = np.zeros(shape=[N_m] * 3)

        # Loop over all mass pairs.  TODO Use bounds or centers here?
        indices, masses = mg.indices(), mg.grid_cell_boundaries()[:-1]
        for i, m_i in zip(indices, masses):
            for j, m_j in zip(indices, masses):
                th = heaviside_theta(i - j)

                # Get collision rate for mass pair.
                R = R_coll[i, j]

                # Calculate combined mass after hit-and-stick collision.
                m_k = m_i + m_j

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
                        K[k_l, i, j] -= R * th * eps
                        K[k_l, i, j] -= R/2 if i==j else 0
                    # Handle "trivial" (non-cancelling) case.
                    else:
                        K[i, i, j] -= R if i<N_m-2 else 0

                # Add "gain" term to kernel.
                # ─────────────────────────────────────────────────────────────
                # Handle upper mass grid boundary.
                k_l = N_m-2 if k_l > N_m-2 else k_l
                k_h = N_m-2 if k_h > N_m-2 else k_h
                if not near_upper_bound:
                    # Handle cancellation.
                    if handle_cancellation:
                        K[k_h, i, j] += R * th * eps
                    # Handle "trivial" (non-cancelling) case.
                    else:
                        K[k_l, i, j] += R * th * (1 - eps)
                        K[k_h, i, j] += R * th * eps

        return K

    def K_frag(self, mg, R_coll):
        indices, masses = mg.indices(), mg.grid_cell_boundaries()
        #                                 ^ TODO Use bounds or centers?

        fragmentation_variants = self.cfg.enable_fragmentation_variant

        K = np.zeros(shape=[mg.N_x] * 3)

        for i, m_i in zip(indices, masses):
            for j, m_j in zip(indices, masses):
                th = heaviside_theta(i - j)

                # 1. "Pulverization"
                # ═════════════════════════════════════════════════════════════════

                if "naive/pulverization" in fragmentation_variants:

                    X = 10
                    k = 0
                    m_k = masses[k]
                    if min(i, j) > X:
                        eps = (m_i + m_j) / m_k
                        K[i, i, j] -= R_coll[i, j]
                        K[k, i, j] += R_coll[i, j] * th * eps

                # todo Implement other variants of fragmentation.
                # - Cratering
                # - Erosion (Watch out for cancellation.)
                # - ...

                # 2. MRN
                # ═════════════════════════════════════════════════════════════════

                if "mrn" in fragmentation_variants:

                    q = -11/6
    
                    # # This is the lowest bin at which fragmentation can occur.
                    X = 10  # (chosen arbitrarily at the moment)
                    if min(i, j) <= X: # Skip fragmentation for masses lower than `X`.
                        continue

                    # Define total mass that needs to be "moved".
                    m_tot = m_i + m_j

                    # Define range of masses resulting from fragmentation event.
                    m_min = masses[0]
                    m_max = m_tot
                    # m_max = max(m_i, m_j)
                    # ^ Note: This is not the maximum allowed mass, but the "first non-allowed".
                    # m_max = max(masses[max(0, i-1)], masses[max(0, j-1)])
                    # ^ Note: Maybe it is necessary to define it like this? (change ">=" to ">" then)

                    A = (q+2) * m_tot / (m_max**(q+2) - m_min**(q+2))  # TODO Calculate from numerical integral

                    # Remove mass from bins corresponding to initial masses.
                    K[i, i, j] -= R_coll[i, j]
                    # ^ Note: Move this into the for loop? (with fraction `f` as prefactor?)

                    # Add mass to bins corresponding to resulting masses.
                    # Loop over all bins that are "receiving" mass.
                    for k, m_k in zip(indices, masses):
                        if m_k >= m_max: # <- Note: See note above. ">=" instead of ">" here.
                            continue

                        # Define fraction `f` of mass `m_tot` landing in bin `k`.
                        n = A * m_k**q
                        f = n / m_tot

                        # Add mass to bin.
                        eps = (m_i + m_j) / m_k
                        K[k, i, j] += R_coll[i, j] * f * th * eps

                # 3. "Splitting"
                # ═════════════════════════════════════════════════════════════════

                # if min(i, j) > X:
                #     l, h = min(i, j), max(i, j)

                #     m_ll = masses[l] / 2
                #     k_l = mass_grid.index_from_value(m_ll)
                #     m_hh = masses[h] / 2
                #     k_h = mass_grid.index_from_value(m_hh)

                #     eps_l = masses[l] / m_ll
                #     eps_h = masses[h] / m_hh

                #     K[l, i, j] -= R[i, j]
                #     K[k_l, i, j] += R[i, j] * th * eps_l
                #     K[h, i, j] -= R[i, j]
                # K[k_h, i, j] += R[i, j] * th * eps_h

        return K
