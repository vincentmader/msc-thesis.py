import numpy as np
from tqdm import tqdm

from axis import DiscreteTimeAxis
from solver.kees_solvers.coag.react_0d import solve_react_0d_equation


SOLVERS = ["explicit_euler", "implicit_euler", "implicit_radau"]


class Solver:

    def __init__(self, cfg):
        self.cfg = cfg
        self.time_axis = DiscreteTimeAxis(cfg)

    def run(self, mg, n_dust, K):
        solver = self.cfg.solver_variant
        N_m = mg.N

        nsubst = 1  # Nr of time substeps between storage of result
        niter = 4   # Nr of iterations for implicit time step

        tg = self.time_axis
        N_t = tg.N
        time = tg.grid_cell_centers  # todo
        # T_START, T_END, N_t = 1e0, 1e9, 200
        # time = np.linspace(0, T_END, N_t)
        # time = T_START * (T_END / T_START)**np.linspace(0, 1, N_t)

        # mgrain = mg.grid_cell_boundaries[:-1]
        # mgrain = m_min * (m_max / m_min)**np.linspace(0, 1, N_m)
        # migrain = np.sqrt(mgrain[1:] * mgrain[:-1])
        # migrain = np.hstack((
        #     migrain[0]**2 / migrain[1],
        #     migrain,
        #     migrain[-1]**2 / migrain[-2]
        # ))
        # dmgrain = migrain[1:] - migrain[:-1]

        # The commented-out lines from above are replaced by calling the methods
        # of the `DiscreteAxis` class instead, this should lead to the same result.
        mgrain = mg.grid_cell_centers
        dmgrain = mg.grid_cell_widths

        # Create initial distribution: Nr of particles per interval dmass per volume.
        # n_dust = np.zeros(N_m)
        # n_dust[0] = 1.0 / dmgrain[0]
        # n_dust[30] = 1.0 / dmgrain[30]
        # Convert this to number of particles per mass bin per volume
        # N_dust = n_dust * dmgrain

        # Convert `n -> N` (number of particles per mass bin per volume).
        N_dust = dmgrain * n_dust

        # def run_solver(N_dust):
        N_dust_store = np.zeros((N_t, N_m))
        N_dust_store[0, :] = N_dust.copy()
        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for itime in tqdm(range(1, N_t)):
            dt = (time[itime] - time[itime - 1]) / nsubst
            for j in range(nsubst):
                assert solver in SOLVERS, f"Unknown solver '{solver}'."
                if solver == "explicit_euler":
                    dNdt = (
                        K[:, :, :] * N_dust[None, :, None] *
                        N_dust[None, None, :]
                    ).sum(axis=2).sum(axis=1)
                    N_dust += dNdt * dt
                if solver == "implicit_euler":
                    N_dust = solve_react_0d_equation(
                        N_dust, s, rmat=rmat, kmat=K,
                        dt=dt, niter=niter, method="backwardeuler"
                    )
                if solver == "implicit_radau":
                    N_dust = solve_react_0d_equation(
                        N_dust, s, rmat=rmat, kmat=K,
                        dt=dt, niter=niter, method="radau"
                    )
                N_dust_store[itime, :] = N_dust.copy()

        # Translate back to physical units
        f = N_dust_store / dmgrain
        m2f = f * mgrain**2
        return N_dust_store, f, m2f


# def dndt(mg, n, K):
#     N_m = mg.N
#     dndt = np.zeros((N_m))
#     for k in range(N_m):
#         dndt_k = 0
#         for i in range(N_m):
#             for j in range(N_m):
#                 dndt_k += K[k][i][j] * n[i] * n[j]
#         dndt[k] = dndt_k
#     return dndt
