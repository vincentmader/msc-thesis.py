import sys

import numpy as np
from tqdm import tqdm

from axis import DiscreteTimeAxis, DiscreteMassAxis
from config import PATH_TO_COAG
from kernel import Kernel, SampledKernel
from visualization.preset.p1 import plot_kernel_gain_loss
from visualization.preset.p2 import plot_kernel_mass_error_vs_time
from visualization.preset.p2 import plot_kernel_sampled_vs_unsampled
from visualization.preset.p2 import plot_sampling_probability_vs_time
from visualization.preset.p2 import plot_sampling_count_vs_time

sys.path.append(PATH_TO_COAG)
from coag.react_0d import solve_react_0d_equation

SOLVERS = ["explicit_euler", "implicit_euler", "implicit_radau"]
N_subst = 1  # Nr of time substeps between storage of result
N_iter  = 4  # Nr of iterations for implicit time step


class Solver:

    def __init__(self, cfg):
        self.cfg = cfg
        self.time_axis = DiscreteTimeAxis(cfg)
        self.mass_axis = DiscreteMassAxis(cfg)

    def run(self, n_dust, K):
        solver = self.cfg.solver_variant

        tg, mg = self.time_axis, self.mass_axis
        tc, mc = tg.bin_centers, mg.bin_centers
        N_t, N_m = tg.N, mg.N
        dm = mg.bin_widths

        # Convert `n -> N` (number of particles per mass bin per volume).
        N_dust = dm * n_dust
        
        R_coag, R_frag = None, None
        kernel = Kernel(self.cfg)
        W_ij = np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)]) 
        #    ^ NOTE Keep definition in sync with W_ij in `SampledKernel`
        R_coll = kernel.R_coag + kernel.R_frag

        Ps, Ns, Ks = [], [], []

        N_dust_store = np.zeros((N_t, N_m))
        N_dust_store[0, :] = N_dust.copy()
        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for itime in tqdm(range(1, N_t)):
            dt = (tc[itime] - tc[itime - 1]) / N_subst

            if self.cfg.enable_collision_sampling:
                if itime == 0:
                    kernel = SampledKernel(self.cfg, N_dust, W_ij=W_ij)
                    R_coag, R_frag = kernel.R_coag, kernel.R_frag
                else:
                    kernel = SampledKernel(self.cfg, N_dust, R_coag=R_coag, R_frag=R_frag, W_ij=W_ij)
                K, P, N = kernel.K, kernel.P_ij, kernel.N_ij
                #    ^ TODO More consistent names, P vs. P_ij
                Ps.append(P)
                Ns.append(N)
                Ks.append(K)

                if False and itime % 20 == 0:
                    plot_kernel_gain_loss(self.cfg, mg, kernel, "log", (1e-20, 1.))

                if False and itime == 160:
                    path = "/Users/vinc/Desktop/K_kij_sampled.txt"
                    np.savetxt(path, kernel.K.reshape(mg.N * mg.N * mg.N))
                    path= "/Users/vinc/Desktop/K_kij_complete.txt"
                    np.savetxt(path, Kernel(self.cfg).K.reshape(mg.N * mg.N * mg.N))
                    path = "/Users/vinc/Desktop/m_k.txt"
                    np.savetxt(path, mc)
                    path = "/Users/vinc/Desktop/dm_k.txt"
                    np.savetxt(path, dm)
                    path = "/Users/vinc/Desktop/n_k.txt"
                    np.savetxt(path, n_dust)

                # Compare kernels: Is sampled with N=2500 same as unsampled?
                if self.cfg.nr_of_samples == self.cfg.mass_resolution**2:
                    kernel_unsampled = Kernel(self.cfg)
                    # Ku = kernel_unsampled.K
                    # err = (K[Ku!=0] - Ku[Ku!=0]) / Ku[Ku!=0]
                    try:
                        # assert (np.abs(err) < 1e-12).all()
                        assert (K == kernel_unsampled.K).all()
                    except AssertionError:
                        # m = (Ku[K!=0] / K[K!=0]).max()
                        # m = ((Ku[K!=0] - K[K!=0]) / K[K!=0]).max()
                        # print(err.max(), err.min())
                        plot_kernel_sampled_vs_unsampled(self.cfg, mg, kernel_unsampled, kernel)

            for _ in range(N_subst):
                assert solver in SOLVERS, f"Unknown solver '{solver}'."
                if self.cfg.solver_variant == "explicit_euler":
                    N_1 = N_dust[None, :, None] # TODO better names than `N_1` and `N_2` ?
                    N_2 = N_dust[None, None, :]
                    dNdt = (K[:, :, :] * N_1 * N_2).sum(axis=2).sum(axis=1)
                    N_dust += dNdt * dt
                if self.cfg.solver_variant == "implicit_euler":
                    N_dust = solve_react_0d_equation(
                        N_dust, s, rmat=rmat, kmat=K,
                        dt=dt, niter=N_iter, method="backwardeuler"
                    )
                if self.cfg.solver_variant == "implicit_radau":
                    N_dust = solve_react_0d_equation(
                        N_dust, s, rmat=rmat, kmat=K,
                        dt=dt, niter=N_iter, method="radau"
                    )
                N_dust_store[itime, :] = N_dust.copy()

        # Translate back to physical units
        f = N_dust_store / dm
        m2f = f * mc**2  # TODO Why multiply with `mgrain`, instead of `dmgrain`?

        dm2f = list(m2f[1:] - m2f[:-1])
        dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
        dm2f = np.array(dm2f)
        dm2f = [dm2f[i] / tg.bin_widths[i] for i, _ in enumerate(dm2f)] # TODO Rename dm2f -> dm2fdt
        # TODO Do the above more elegantly. (Calculate temp. deriv.)

        if self.cfg.enable_collision_sampling:
            plot_sampling_probability_vs_time(self.cfg, mg, Ps)
            plot_sampling_count_vs_time(self.cfg, mg, Ns)
            plot_kernel_mass_error_vs_time(self.cfg, mg, Ks, R_coll)

        return N_dust_store, f, m2f, dm2f
