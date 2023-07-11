import numpy as np
from tqdm import tqdm

from axis import DiscreteTimeAxis
from solver.kees_solvers.coag.react_0d import solve_react_0d_equation


SOLVERS = ["explicit_euler", "implicit_euler", "implicit_radau"]
N_subst = 1  # Nr of time substeps between storage of result
N_iter = 4  # Nr of iterations for implicit time step


class Solver:

    def __init__(self, cfg):
        self.cfg = cfg
        self.time_axis = DiscreteTimeAxis(cfg)

    def run(self, mg, n_dust, K):
        solver = self.cfg.solver_variant
        N_m = mg.N

        N_t = self.time_axis.N
        time = self.time_axis.grid_cell_centers 

        # The commented-out lines from above are replaced by calling the methods
        # of the `DiscreteAxis` class instead, this should lead to the same result.
        mgrain = mg.grid_cell_centers
        dmgrain = mg.grid_cell_widths

        # Convert `n -> N` (number of particles per mass bin per volume).
        N_dust = dmgrain * n_dust

        N_dust_store = np.zeros((N_t, N_m))
        N_dust_store[0, :] = N_dust.copy()
        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for itime in tqdm(range(1, N_t)):
            dt = (time[itime] - time[itime - 1]) / N_subst
            for j in range(N_subst):
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
                        dt=dt, niter=N_iter, method="backwardeuler"
                    )
                if solver == "implicit_radau":
                    N_dust = solve_react_0d_equation(
                        N_dust, s, rmat=rmat, kmat=K,
                        dt=dt, niter=N_iter, method="radau"
                    )
                N_dust_store[itime, :] = N_dust.copy()

        # Translate back to physical units
        f = N_dust_store / dmgrain
        m2f = f * mgrain**2  # TODO Why multiply with `mgrain`, instead of `dmgrain`?
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
