import numpy as np

from config import Config


def collision_outcome_probabilities(cfg: Config, dv: np.ndarray):
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
    return P_coag, P_frag


def collision_outcome_probabilities_from_maxwell_boltzmann(
    cfg: Config,
    dv,
):
    v_frag = cfg.fragmentation_velocity
    N_m = cfg.mass_resolution
    P_coag = np.zeros([N_m]*2)
    P_frag = np.zeros([N_m]*2)

    for i in range(N_m):
        for j in range(N_m):
            # Calculate fragmentation probability.
            P_f = (3/2 * (v_frag/dv[i, j])**2 + 1) * np.exp(-3/2 * (v_frag/dv[i, j])**2)
            # ^ See "2022 Stammler & Birnstiel"
            P_frag[i, j] = P_f 
            P_coag[i, j] = 1 - P_f  # NOTE: Bouncing is neglected here.

    P_tot = np.sum(P_coag + P_frag) / N_m**2
    assert P_tot == 1

    return P_coag, P_frag


def collision_outcome_probabilities_from_cutoff_velocity(
    cfg: Config,
    dv,
):
    v_frag = cfg.fragmentation_velocity
    N_m = cfg.mass_resolution
    P_coag = np.zeros([N_m]*2)
    P_frag = np.zeros([N_m]*2)

    for i in range(N_m):
        for j in range(N_m):
            # Decide whether fragmentation will occur.
            will_fragment = dv[i, j] > v_frag
            if will_fragment:
                P_frag[i, j] = 1
            else:  # NOTE: Bouncing is neglected here.
                P_coag[i, j] = 1

    P_tot = np.sum(P_coag + P_frag) / N_m**2
    assert P_tot == 1

    return P_coag, P_frag
