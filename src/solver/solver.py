from pathlib import Path
import sys

import numpy as np
from tqdm import tqdm

from axis import DiscreteTimeAxis
from config import PATH_TO_LIB
path = Path(PATH_TO_LIB, "coag_py")
path = str(path)
sys.path.append(path)
from coag.react_0d import solve_react_0d_equation
from dust import particle_radius_from_mass
from sampling.kernel import SampledKernel
from visualization.base import GridspecPlot, PcolorMatrixSubplot
from visualization.v1.mass_conservation import KernelMassConservationPlot

SOLVERS = ["explicit_euler", "implicit_euler", "implicit_radau"]
N_subst = 1  # Nr of time substeps between storage of result
N_iter = 4  # Nr of iterations for implicit time step

def plot(kernel): # TODO Move elsewhere.
    cfg, mg, K = kernel.cfg, kernel.mg, kernel.K
    ac = particle_radius_from_mass(mg.grid_cell_centers, cfg.dust_particle_density)
    GridspecPlot([
        PcolorMatrixSubplot(
            ac, ac, K,
            title="kernel gain contribution $G_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            symmetrized=True,
            z_limits=(1e-20, 1e-7),
        ),
        PcolorMatrixSubplot(
            ac, ac, K,
            title="kernel gain contribution $G_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            symmetrized=True,
            z_limits=(1e-20, 1e-7),
        )
    ], add_slider=True).render()

    p = KernelMassConservationPlot(cfg, mg, K)
    p.show()


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

            if self.cfg.enable_collision_sampling:
                kernel = SampledKernel(self.cfg, N_dust)
                K = kernel.K
                if itime % 10 == -1:
                    plot(kernel) 

            for j in range(N_subst):
                assert solver in SOLVERS, f"Unknown solver '{solver}'."
                if solver == "explicit_euler":
                    N_1 = N_dust[None, :, None] # TODO better names than `N_1` and `N_2` ?
                    N_2 = N_dust[None, None, :]
                    dNdt = (K[:, :, :] * N_1 * N_2).sum(axis=2).sum(axis=1)
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
