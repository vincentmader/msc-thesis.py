from dataclasses import field
import sys

import numpy as np
from tqdm import tqdm

from config import Config
from config import PATH_TO_COAG, SOLVERS
from functions.disk import mass_distribution
from functions.plotting.preset.p2 import plot_kernel_mass_error_vs_time
from functions.plotting.preset.p2 import plot_kernel_sampled_vs_unsampled
from functions.plotting.preset.p2 import plot_sampling_probability_vs_time
from functions.plotting.preset.p2 import plot_sampling_count_vs_time
from functions.solver import fix_negative_densities
from models.axis import DiscreteMassAxis, DiscreteRadialAxis, DiscreteTimeAxis
from models.disk import Disk, DiskRegion
from models.kernel import Kernel, SampledKernel

sys.path.append(PATH_TO_COAG)
from coag.react_0d import solve_react_0d_equation


N_subst = 1  # Nr of time substeps between storage of result
N_iter  = 4  # Nr of iterations for implicit time step


class Solver():

    cfg:                    Config
    disk:                   Disk
    region:                 DiskRegion
    rg:                     DiscreteRadialAxis
    mg:                     DiscreteMassAxis
    tg:                     DiscreteTimeAxis

    N_k_vs_t:               np.ndarray 
    n_k_vs_t:               np.ndarray
    M_k_vs_t:               np.ndarray
    dMdt_vs_t:              np.ndarray

    P_ij_vs_t:              np.ndarray
    S_ij_vs_t:              np.ndarray
    dK_ij_vs_t:             np.ndarray

    kernels:                list[Kernel]    = field(default_factory=list)
    iteration_idx:          int             = 0

    def __init__(self, cfg: Config):

        self.cfg        = cfg
        self.disk       = Disk(cfg)
        self.region     = DiskRegion(cfg, self.disk)
        self.rg         = DiscreteRadialAxis(cfg)
        self.mg         = DiscreteMassAxis(cfg)
        self.tg         = DiscreteTimeAxis(cfg)

        self.n_k_vs_t   = np.zeros(shape=(self.tg.N, self.mg.N))
        self.N_k_vs_t   = np.zeros(shape=(self.tg.N, self.mg.N))
        self.M_k_vs_t   = np.zeros(shape=(self.tg.N, self.mg.N))
        self.dMdt_vs_t  = np.zeros(shape=(self.tg.N, self.mg.N))

        self.P_ij_vs_t  = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))
        self.S_ij_vs_t  = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))
        self.dK_ij_vs_t = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))

        self.kernels    = [Kernel(cfg)]

    def run(self):

        N_m, N_t    = self.mg.N,            self.tg.N
        mc,  tc     = self.mg.bin_centers,  self.tg.bin_centers
        dm          = self.mg.bin_widths

        self.N_k_vs_t[0, :] = dm * mass_distribution.dirac_delta(cfg)

        W_ij = np.sum([
            mc[k] * np.abs(self.kernels[0].K[k]) for k in range(N_m)
        ])

        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for i_t in tqdm(range(1, N_t)):
            self.iteration_idx = i_t

            if self.cfg.enable_collision_sampling:
                sampled = SampledKernel(cfg, self.N_k_vs_t[i_t-1], W_ij=W_ij)
                self.kernels.append(sampled)
                self.P_ij_vs_t[i_t] = sampled.P_ij
                self.S_ij_vs_t[i_t] = sampled.N_ij  # TODO Rename.

            self.step(rmat, s)

        n    = self.N_k_vs_t * dm
        M    = self.N_k_vs_t * mc
        dMdt = np.array([
            (M[i] - M[i-1]) / (tc[i] - tc[i-1]) for i in range(1, self.tg.N)
        ])
        self.n_k_vs_t       = n
        self.M_k_vs_t       = M
        self.dMdt_vs_t[1:]  = dMdt

    def step(self, rmat, s):
        i_t = self.iteration_idx
        tc = self.tg.bin_centers
        dt = (tc[i_t] - tc[i_t - 1]) / N_subst
        K = self.kernels[-1].K

        N_dust = self.N_k_vs_t[i_t - 1]

        for _ in range(N_subst):
            assert self.cfg.solver_variant in SOLVERS, f"{self.cfg.solver_variant} = "

            if self.cfg.solver_variant == "explicit_euler":
                N_i = N_dust[None, :, None]
                N_j = N_dust[None, None, :]
                dNdt = (K[:, :, :] * N_i * N_j).sum(axis=2).sum(axis=1)
                N_dust += dNdt * dt

            if self.cfg.solver_variant == "implicit_euler":
                N_dust = solve_react_0d_equation(N_dust, s, 
                    rmat=rmat, kmat=K, dt=dt, niter=N_iter, method="backwardeuler"
                )

            if self.cfg.solver_variant == "implicit_radau":
                N_dust = solve_react_0d_equation(N_dust, s, 
                    rmat=rmat, kmat=K, dt=dt, niter=N_iter, method="radau"
                )

        N_dust = fix_negative_densities(self.mg.bin_centers, N_dust)
        self.N_k_vs_t[i_t, :] = N_dust.copy()


if __name__ == "__main__":

    cfg     = Config(
        enable_collision_sampling=True,
        # nr_of_samples=150,
    )
    solver = Solver(cfg)
    solver.run()
    
    N       = solver.N_k_vs_t
    n       = solver.n_k_vs_t
    M       = solver.M_k_vs_t
    dMdt    = solver.dMdt_vs_t
    
    tc      = solver.tg.bin_centers
    mg      = solver.mg
    dm      = mg.bin_widths
    
    kernel  = solver.kernels[-1]  # TODO
    
    from functions.plotting.preset.p1 import plot_evolution
    plot_evolution(cfg, mg, "", tc, N, n, M, dMdt)
    from functions.plotting.preset.p1 import plot_error
    plot_error(cfg, mg, kernel, tc, np.sum(M, axis=1))

    if cfg.enable_collision_sampling:

        P_ij_vs_t = solver.P_ij_vs_t
        plot_sampling_probability_vs_time(cfg, mg, P_ij_vs_t, symmetrized=False)

        S_ij_vs_t = solver.S_ij_vs_t
        plot_sampling_count_vs_time(cfg, mg, S_ij_vs_t, symmetrized=False)

        kernel = solver.kernels[0]
        R_coll = kernel.R_coag + kernel.R_frag
        K_kij_vs_t = [k.K for k in solver.kernels]
        plot_kernel_mass_error_vs_time(cfg, mg, K_kij_vs_t, R_coll)

        # Ns = [kernel.N_ij for kernel in self.kernels]
        print("\"Actual\" sampling density:")
        print("N_sample_tot / (N_t * N_m^2) * 2 =", np.sum(S_ij_vs_t) / solver.tg.N / mg.N**2 * 100 * 2, "%")

        # TODO Calculate: 
        #      - nr. of kernel entries == 0
        #      - nr. of kernel entries << 1
        #      - nr. of kernel entriee ~~ 1
        #      - nr. of kernel entries sampled

#     # Compare kernels: Is sampled with N=2500 same as unsampled?
#     if self.cfg.nr_of_samples == self.cfg.mass_resolution**2:
#         kernel_unsampled = Kernel(self.cfg)
#         try:
#             assert (K == kernel_unsampled.K).all()
#         except AssertionError:
#             plot_kernel_sampled_vs_unsampled(self.cfg, self.mg, kernel_unsampled, kernel)
