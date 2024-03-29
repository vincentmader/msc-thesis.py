import sys

import numpy as np
from tqdm import tqdm

from axis import DiscreteTimeAxis, DiscreteMassAxis, AxisLabelVariant
from config import PATH_TO_COAG
from kernel import Kernel, SampledKernel
from visualization.base import GridspecPlot
from visualization.kernel import KernelSubplot, KernelMassConservationSubplot
from visualization.preset.p1 import plot_kernel_gain_loss

sys.path.append(PATH_TO_COAG)
from coag.react_0d import solve_react_0d_equation

SOLVERS = ["explicit_euler", "implicit_euler", "implicit_radau"]
N_subst = 1  # Nr of time substeps between storage of result
N_iter = 4  # Nr of iterations for implicit time step

def plot(cfg, mg, K, P): # TODO Move elsewhere.
    GridspecPlot([
        KernelSubplot(
            cfg, mg, K, axis_label_variant=AxisLabelVariant.Radius,
            title="kernel gain contribution $G_{kij}$",
            symmetrized=True,
            z_limits=(1e-20, 1e-7),
        ),
        KernelSubplot(
            cfg, mg, K, axis_label_variant=AxisLabelVariant.Radius,
            title="kernel gain contribution $G_{kij}$",
            symmetrized=True,
            z_limits=(1e-20, 1e-7),
        ),
    ], add_slider=True).render()

    GridspecPlot([
        KernelSubplot(
            cfg, mg, P, title="sampling probability $P_{ij}$",
        )
    ]).render()

    GridspecPlot([
        KernelMassConservationSubplot(cfg, mg, K)
    ]).render()


def plot_2(cfg, mg, Ps):
    GridspecPlot(
        [
            KernelSubplot(
                cfg, mg, np.array(Ps), cmap="Blues", z_limits=(1e-5, 1),
                title="Collision Pair Sampling Probability $P_{ij}$",
            ),
        ], add_slider=True,
    ).render()


def plot_3(cfg, mg, Ns):
    GridspecPlot(
        [
            KernelSubplot(
                cfg, mg, np.array(Ns), cmap="Blues", 
                z_limits=(0, np.max(Ns)),  # <- NOTE: Low upper boundary for better visibility.
                axis_scales=("log", "log", "lin"),
                title=r"Collision Pair Sampling Count $N_{ij}$",
            ),
        ], add_slider=True,
    ).render()


def plot_4(cfg, mg, kernel_unsampled, kernel_sampled):
    K_g = kernel_unsampled.K_gain
    K_l = kernel_unsampled.K_loss
    S_g = kernel_sampled.K_gain
    S_l = kernel_sampled.K_loss

    cmap = "Reds" 
    z_limits = (1e-20, 1)
    scale = "log"
    GridspecPlot([
        KernelSubplot(
            cfg, mg, K_g,
            title="kernel gain contribution $G_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            symmetrized=False,
            cmap=cmap,
        ),
        KernelSubplot(
            cfg, mg, -K_l,
            title="kernel loss contribution $L_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            ylabel="",
            symmetrized=False,
            cmap=cmap,
        ),
        KernelSubplot(
            cfg, mg, S_g,
            title="kernel gain contribution $G_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            symmetrized=False,
            cmap=cmap,
        ),
        KernelSubplot(
            cfg, mg, -S_l,
            title="kernel loss contribution $L_{kij}$",
            axis_scales=(scale, scale, scale),
            z_limits=z_limits,
            ylabel="",
            symmetrized=False,
            cmap=cmap,
        ),
    ], add_slider=True).render()



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
        W_ij = np.sum([mc[k] * np.abs(K[k]) for k in range(mg.N)]) # NOTE Keep in sync with W_ij in `SampledKernel`

        Ps, Ns = [], []

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

                kernel_unsampled = Kernel(self.cfg)
                if self.cfg.nr_of_samples == self.cfg.mass_resolution**2:
                    try:
                        assert (K == kernel_unsampled.K).all()
                    except AssertionError:
                        plot_4(self.cfg, mg, kernel_unsampled, kernel)
                # TODO Compare kernels: Is sampled with N=2500 same as unsampled?

                # K = [(K_k + K_k.T)/2 for K_k in K]  # TODO Does this make difference for integration?

                # for k, K_k in enumerate(K):  # NOTE This cannot be correct! Divide by N_ij instead!
                #     K[k][P!=0] = K[k][P!=0] / P[P!=0]

                # for k, K_k in enumerate(K):  
                #     K[k][N_ij != 0] = K[k][N_ij != 0] / N_ij[N_ij != 0]
                # for k, K_k in enumerate(K):  
                #     assert np.sum(N_ij) == self.cfg.nr_of_samples
                #     K[k] = K[k] / np.sum(N_ij)

                # for i in range(self.cfg.mass_resolution):
                #     for j in range(self.cfg.mass_resolution):
                #         S = sum([K[k,i,j] for k in range(self.cfg.mass_resolution)])
                #         if S != 0:
                #             print(i, j, S)

                # K = [.5 * (K_k + K_k.T) for K_k in K]
                # if itime % 10 == 0:
                # if itime in [1, 50, 100, 120, 130, 140]:
                # if itime > 150 and itime % 5 == 0:
                #     plot(self.cfg, mg, K, P)
                #     a = input()
                #     if a == "p":
                #         path = Path(f"/Users/vinc/Desktop/K_kij_{itime}.txt")
                #         kernel.save_to_file(path=path)

            for _ in range(N_subst):
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
        f = N_dust_store / dm
        m2f = f * mc**2  # TODO Why multiply with `mgrain`, instead of `dmgrain`?

        dm2f = list(m2f[1:] - m2f[:-1])
        dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
        dm2f = np.array(dm2f)
        dm2f = [dm2f[i] / tg.bin_widths[i] for i, _ in enumerate(dm2f)] # TODO Rename dm2f -> dm2fdt
        # TODO Do the above more elegantly. (Calculate temp. deriv.)

        if self.cfg.enable_collision_sampling:
            plot_2(self.cfg, mg, Ps)
            plot_3(self.cfg, mg, Ns)

        return N_dust_store, f, m2f, dm2f
