from dataclasses import dataclass
import sys

import numpy as np
from tqdm import tqdm

from models.axis import DiscreteTimeAxis, DiscreteMassAxis
from config import Config
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


@dataclass
class Solver:
    cfg:            Config
    time_axis:      DiscreteTimeAxis
    mass_axis:      DiscreteMassAxis

    # N_dust:         np.ndarray
    # N_dust_vs_t:    np.ndarray

    kernels:        list[SampledKernel] | list[Kernel] | list[np.ndarray]

    def __init__(
        self, 
        cfg: Config,
    ):

        self.cfg        = cfg
        self.time_axis  = DiscreteTimeAxis(cfg)
        self.mass_axis  = DiscreteMassAxis(cfg)
        self.kernels    = []

    def run(self, n_dust, K):

        tg, mg = self.time_axis, self.mass_axis
        tc, mc = tg.bin_centers, mg.bin_centers
        N_t, N_m = tg.N, mg.N
        dm = mg.bin_widths

        # Convert `n -> N` (number of particles per mass bin per volume).
        N_dust = dm * n_dust
        N_dust_store = np.zeros((N_t, N_m))
        N_dust_store[0, :] = N_dust.copy()
        
        kernel = Kernel(self.cfg)
        W_ij = np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)]) 
        #    ^ NOTE Keep definition in sync with W_ij in `SampledKernel`
        R_coll = kernel.R_coag + kernel.R_frag

        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for itime in tqdm(range(1, N_t)):
            dt = (tc[itime] - tc[itime - 1]) / N_subst

            if self.cfg.enable_collision_sampling:
                kernel = SampledKernel(self.cfg, N_dust, W_ij=W_ij)
                self.kernels.append(kernel)
                K = kernel.K
                # K = [(K_k + K_k.T)/2 for K_k in K]
                # TODO Find out: Do I HAVE to be symmetrize?

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
                    try:
                        assert (K == kernel_unsampled.K).all()
                    except AssertionError:
                        plot_kernel_sampled_vs_unsampled(self.cfg, mg, kernel_unsampled, kernel)

            for _ in range(N_subst):
                assert self.cfg.solver_variant in SOLVERS, f"{self.cfg.solver_variant} = "
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

            N_dust = fix_negative_densities(mc, N_dust)

        # Translate back to physical units
        f   = N_dust_store / dm
        m2f = N_dust_store * mc
        # m2f = f * mc**2  # NOTE Why multiply with `mgrain`, instead of `dmgrain`?

        # Calculate temporal derivative.
        dm2f_0 = np.zeros(shape=(mg.N))
        dm2f = [dm2f_0] + [m2f[i] - m2f[i-1] for i in range(1, tg.N)]
        dm2fdt  = [dm2f[i] / tg.bin_widths[i] for i, _ in enumerate(tc)]

        if self.cfg.enable_collision_sampling:
            if False:
                Ps = [kernel.P_ij for kernel in self.kernels]
                plot_sampling_probability_vs_time(self.cfg, mg, Ps, symmetrized=False)

            if True:
                Ns = [kernel.N_ij for kernel in self.kernels]
                print("\"Actual\" sampling density:")
                print("N_sample_tot / (N_t * N_m^2) * 2 =", np.sum(Ns) / tg.N / mg.N**2 * 100 * 2, "%")
                # TODO Calculate: 
                #      - nr. of kernel entries == 0
                #      - nr. of kernel entries << 1
                #      - nr. of kernel entriee ~~ 1
                #      - nr. of kernel entries sampled

            if False:
                plot_sampling_count_vs_time(self.cfg, mg, Ns)

            if False:
                Ks = [kernel.K for kernel in self.kernels]
                plot_kernel_mass_error_vs_time(self.cfg, mg, Ks, R_coll)

        return N_dust_store, f, m2f, dm2fdt

def fix_negative_densities(
    mc:     np.ndarray,
    N_dust: np.ndarray,
):
    idx_i = np.where(N_dust < 0)
    idx_k = np.where(N_dust >= 0)

    M_i = np.sum(N_dust[idx_i] * mc[idx_i])
    M_k = np.sum(N_dust[idx_k] * mc[idx_k])

    fac = np.abs(M_i) / np.abs(M_k)

    N_dust[idx_k] *= 1 - fac
    N_dust[idx_i] = 0

    return N_dust

def deriv(
    y: np.ndarray,
    x: np.ndarray,
):
    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    return dy / dx
