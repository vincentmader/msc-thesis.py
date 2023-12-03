import numpy as np

from config import Config
from functions.plotting.preset.p1 import plot_error
from functions.plotting.preset.p1 import plot_evolution
from functions.plotting.preset.p2 import plot_kernel_mass_error_vs_time
from functions.plotting.preset.p2 import plot_sampling_count_vs_time
from functions.plotting.preset.p2 import plot_sampling_probability_vs_time
from models.solver import SolverV2

if __name__ == "__main__":

    cfg     = Config(enable_collision_sampling=True)
    solver = SolverV2(cfg)
    solver.run()
    
    N       = solver.N_k_vs_t
    n       = solver.n_k_vs_t
    M       = solver.M_k_vs_t
    dMdt    = solver.dMdt_vs_t
    tc      = solver.tg.bin_centers
    mg      = solver.mg
    dm      = mg.bin_widths
    
    plot_evolution(cfg, mg, "", tc, N, n, M, dMdt)
    plot_error(cfg, mg, tc, np.sum(M, axis=1))

    if cfg.enable_collision_sampling:

        K_kij_vs_t = [k.K for k in solver.kernels]
        P_ij_vs_t = solver.P_ij_vs_t
        S_ij_vs_t = solver.S_ij_vs_t

        plot_sampling_probability_vs_time(cfg, mg, P_ij_vs_t, symmetrized=False)
        plot_sampling_count_vs_time(cfg, mg, S_ij_vs_t, symmetrized=False)
        plot_kernel_mass_error_vs_time(cfg, mg, K_kij_vs_t, solver.kernels[0].R_coll)

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
