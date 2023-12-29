import sys

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from config import Config
from config import PATH_TO_COAG
from functions.plotting.preset.p1 import plot_error
from functions.plotting.preset.p1 import plot_evolution
from functions.plotting.preset.p2 import plot_kernel_mass_error_vs_time
from functions.plotting.preset.p2 import plot_sampling_count_vs_time
from functions.plotting.preset.p2 import plot_sampling_probability_vs_time
from models.kernel import Kernel
from models.kernel import SampledKernel
# from models.kernel import SampledKernelV2
from models.solver import SolverV2


if __name__ == "__main__":

    kernel = Kernel(Config(enable_collision_sampling=False)) 
    cfg    = Config(enable_collision_sampling=False)  # note: Only change this, not above.
    solver = SolverV2(cfg)
    solver.run()
    
    N       = solver.N_k_vs_t
    n       = solver.n_k_vs_t
    M       = solver.M_k_vs_t
    dMdt    = solver.dMdt_vs_t
    tg      = solver.tg
    tc      = tg.bin_centers
    mg      = solver.mg
    mc      = mg.bin_centers
    dm      = mg.bin_widths
    
    plot_evolution(cfg, mg, "", tc, N, n, M, dMdt)
    # plot_error(cfg, mg, tc, np.sum(M, axis=1))

    if cfg.enable_collision_sampling:

        K_kij_vs_t = [k.K for k in solver.kernels]
        P_ij_vs_t = solver.P_ij_vs_t
        S_ij_vs_t = solver.S_ij_vs_t

        # plot_sampling_probability_vs_time(cfg, mg, P_ij_vs_t, symmetrized=False)
        # plot_sampling_count_vs_time(cfg, mg, S_ij_vs_t, symmetrized=False)
        # plot_kernel_mass_error_vs_time(cfg, mg, K_kij_vs_t, solver.kernels[0].R_coll)

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

    sys.path.append(PATH_TO_COAG)
    from coag.react_0d import solve_react_0d_equation
    N_iter  = 4        # = nr. of iterations for implicit time step
    N_subst = 1        # = nr. of time substeps between storage of result
    s = np.zeros((mg.N))
    rmat = np.zeros((mg.N, mg.N))

    N_avg    = 150      # = nr. of samples to average derivative over
    i_t      = 180     # = iteration index (time) 
    N_t      = N[i_t]  # = mass distribution at that time
    dt = tc[i_t] - tc[i_t-1]  # = time step size at that time

    # Pre-calculate weights for later definition of sampling probability.
    cfg    = Config(enable_collision_sampling=True)
    W_ij = SampledKernel(cfg, N_t, N_sample=mg.N**2).W_ij
    # W_ij = SampledKernelV2(cfg, N_t, N_sample=mg.N**2).W_ij
    # assert W_ij_1.all() == W_ij.all()

    # Using complete kernel, forward the mass distribution one step.
    N_complete = solve_react_0d_equation(N_t.copy(), s, 
        rmat=rmat, kmat=kernel.K, dt=dt/N_subst, niter=N_iter, method="backwardeuler"
    )
    # Calculate total change, & num. derivative.
    dN_complete = N_complete - N_t
    dNdt_complete = dN_complete / dt

    # For each sampling density value, define nr. of to-sample kernel entries.
    # sampling_densities = [ 0.2, 0.4, 0.6, 0.8, 1.0 ]
    sampling_densities = [ 0.6, 0.8, 1.0 ]
    for i, sampling_density in enumerate(sampling_densities):
        N_sample = int((sampling_density * (mg.N**2 - mg.N) / 2 + mg.N)) # TODO Fix parentheses
        # todo Define sampling density with or without 1/2 ?

        dN = []
        for batch_idx in tqdm(range(N_avg)):

            # Build sampled kernel.
            kernel = SampledKernel( cfg, N_t, N_sample=N_sample, W_ij=W_ij )
            K = kernel.K

            # From kernel, calculate change in mass distribution.
            N_tp1 = solve_react_0d_equation(N_t.copy(), s, 
                rmat=rmat, kmat=K, dt=dt/N_subst, niter=N_iter, method="backwardeuler"
            )
            assert np.sum( (N_tp1 - N_t) * mc * dm, axis=0 ) < 1e-16
            dN.append(N_tp1 - N_t)  # Note: `tp1` = `t + 1`.

        # Plot.
        NR_OF_BATCHES_DISPLAYED = 4
        COLORS = ["r", "g", "b", "k", "c", "m", "y"]
        plot_idx = 0
        plot_at = [int(i) for i in np.linspace(0, 1, NR_OF_BATCHES_DISPLAYED) * N_avg] + [NR_OF_BATCHES_DISPLAYED-1]
        for batch_idx in range(N_avg):
            if batch_idx not in plot_at:
                continue

            dN_avg = sum([dN for dN in dN[:batch_idx+1]]) / (batch_idx + 1)
            # y = (dN_avg*mc*dm - dN_complete*mc*dm) / (dN_complete*mc*dm)
            y = (dN_avg - dN_complete) / dN_complete
            plt.semilogy(y,  f"{COLORS[plot_idx]}-", label=r"$\rho_s=$" + f"{sampling_density}, sampled {batch_idx+1} times")
            plt.semilogy(-y, f"{COLORS[plot_idx]}--")
            plot_idx += 1

        plt.xlabel("mass bin index")
        ylabel = r"rel. diff.: $\frac{\dot M_{sampled} - \dot M_{complete}}{\dot M_{complete}}$"
        plt.ylabel(ylabel)
        plt.legend()
        plt.show()
        plt.close()
