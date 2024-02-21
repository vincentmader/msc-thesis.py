import os, sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
sys.path.append(os.path.join("..", "..", "src"))
from config import Config
from config import Config
from config import PATH_TO_FIGURES, PATH_TO_OUTFILES
from config import PATH_TO_COAG
from constants import AU
from functions.dust.collision import collision_outcome_probabilities_from_maxwell_boltzmann
from functions.dust.collision import collision_rate
from functions.dust.relative_velocity import relative_velocity
from functions.utils.dates import format_seconds_as_years
from functions.solver import fix_negative_densities
from models.axis import DiscreteMassAxis, DiscreteRadialAxis
from models.axis import DiscreteTimeAxis
from models.disk import Disk, DiskRegion
from models.kernel import Kernel
from models.kernel import SampledKernel
from models.solver import SolverV2
sys.path.append(PATH_TO_COAG)
from coag.react_0d import solve_react_0d_equation

path_to_outfiles = Path(PATH_TO_OUTFILES, "data", "113")
path_to_figures  = Path(PATH_TO_FIGURES, "113")
os.makedirs(path_to_figures, exist_ok=True)


# Define model variables.
# ═════════════════════════════════════════════════════════════════════════════
N_m = 50
i_t = 198

rho_sample = 0.5
nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)

N_subst = 1
N_iter = 4
rmat = np.zeros((N_m, N_m))
s = np.zeros((N_m))

# Integrate to equilibrium using full kernel.
# ═════════════════════════════════════════════════════════════════════════════

cfg = Config(
    mass_resolution=N_m,
    enable_collision_sampling=False,
)
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
tg = DiscreteTimeAxis(cfg)
tc = tg.bin_centers
mc = mg.bin_centers
dt = tc[i_t+1] - tc[i_t]

disk = Disk(cfg, rg, mg)
disk_region = DiskRegion(cfg, disk)
dv = relative_velocity(cfg, disk, disk_region)
R_coll = collision_rate(cfg, disk, disk_region)
P_coag, P_frag = collision_outcome_probabilities_from_maxwell_boltzmann(cfg, dv)
R_coag = R_coll * P_coag
R_frag = R_coll * P_frag

solver = SolverV2(cfg)
solver.run()

N_eq    = solver.N_k_vs_t[i_t]
M_eq    = solver.M_k_vs_t[i_t]
dM_eq   = solver.M_k_vs_t[i_t+1] - solver.M_k_vs_t[i_t]
dMdt_eq = dM_eq / dt

# Sample the kernel from complete solution a couple of times
# ═════════════════════════════════════════════════════════════════════════════

cfg = Config(
    mass_resolution=N_m,
    nr_of_samples=nr_of_samples,
    enable_collision_sampling=True,
)

N_s = 30
M_vs_s = []
for i_s in tqdm(range(N_s)):
    kernel = SampledKernel(cfg, N_eq, R_coag, R_frag)

    N = solve_react_0d_equation(N_eq, s, 
        rmat=rmat, kmat=kernel.K, dt=dt, niter=N_iter, method="radau"
    )
    N = fix_negative_densities(mg.bin_centers, N)
    M = N * mc
    M_vs_s.append(M)


y = M_vs_s
y = [(y[i_s] - M_eq) / M_eq for i_s in range(N_s)]
y = [(y[i_s]**2).sum()**.5 for i_s in range(N_s)]
# y = M_vs_s
# y = [y[i_s] - M_eq for i_s in range(N_s)]  # dM
# y = [y[i_s] / dt for i_s in range(N_s)]    # dMdt
# y = [sum(y[:i_s])/i_s for i_s in range(1, N_s)]
# y = [(y[i_s] - dMdt_eq) / dMdt_eq for i_s in range(N_s-1)]
# y = [(y[i_s]**2).sum()**.5 for i_s in range(N_s-1)]
# X = (dMdt_eq**2).sum()**.5
# y = [(y[i_s] - X) / X for i_s in range(N_s-1)]

plt.semilogy(y)
plt.show()
plt.close()


# y = 

# Plot
# ═════════════════════════════════════════════════════════════════════════════

# import numpy as np
# from matplotlib.widgets import Slider

# def update(val):
#     i = int(val)
#     line1.set_ydata(M_vs_s[i]) 
#     line2.set_ydata(dMdt_vs_s[i]) 

# plt.subplot(2, 1, 1)
# line1, = plt.semilogy(M_vs_s[0]) 
# plt.subplot(2, 1, 2)
# line2, = plt.semilogy(dMdt_vs_s[0]) 

# # plt.ylim(1e-25, 1)

# freq_slider_ax = plt.axes([0.25, 0.1, 0.65, 0.03])  # Define position and size of the slider
# freq_slider = Slider(freq_slider_ax, '', 0.1, 10.0, valinit=0)  # Create the slider
# freq_slider.on_changed(update)  # Call update function when slider value changes

# plt.show()
# plt.close()

# # print(dMdt_vs_s)
# # print(len(dMdt_vs_s))
# # y = [dMdt_vs_s[i_s] for i_s in range(1, N_s)]
# for i_s in range(N_s):
#     plt.plot(dMdt_vs_s[i_s], label=i_s)
# plt.legend()
# # y = [i.sum() for i in y]
# # plt.plot(y)
# plt.show()
# plt.close()





# # N_dust = solver.N_k_vs_t[i_t]
# # dN_dust = (solver.N_k_vs_t[i_t+1] - solver.N_k_vs_t[i_t]) / (tc[i_t+1] - tc[i_t])

# # cfg.enable_collision_sampling = True

# # N_dusts = [N_dust]
# # for i in tqdm(range(N_batch)):

# #     kernel = SampledKernel(cfg, N_dust, R_coag, R_frag)
    
# #     dt = (tc[i_t] - tc[i_t - 1]) / N_subst
    
# #     N = solve_react_0d_equation(N_dust, s, 
# #         rmat=rmat, kmat=kernel.K, dt=dt, niter=N_iter, method="radau"
# #     )
# #     N = fix_negative_densities(mg.bin_centers, N)
# #     N_dusts.append(N)

# # y = [N_dusts[i] - N_dust for i in range(1, N_batch)]
# # # y = [i / (tc)]
# # y = [((i**2).sum())**.5 for i in y]
# # y = [sum(y[:i]) / i for i in range(1, len(y)) ]
# # plt.semilogy(y)
# # plt.show()
# # plt.close()


# # i_t = 20

# # plt.subplot(3, 1, 1)
# # plt.loglog(mc, N_complete)
# # plt.subplot(3, 1, 2)
# # plt.loglog(mc, dN_complete)
# # plt.subplot(3, 1, 3)
# # plt.loglog(mc, dNdt_complete)
# # plt.show()
# # plt.close()

# # N_s = 20
# # dNdts_sampled = []
# # for i_s in tqdm(range(N_s)):
# #     kernel = SampledKernel(cfg, N_complete, R_coag, R_frag)

# #     N_sampled = solve_react_0d_equation(N_complete, s, 
# #         rmat=rmat, kmat=kernel.K, dt=dt, niter=N_iter, method="radau"
# #     )
# #     N_sampled = fix_negative_densities(mg.bin_centers, N_sampled)

# #     dN_sampled = N_sampled - N_complete
# #     dNdt_sampled = dN_sampled / dt

# #     dNdts_sampled.append(dNdt_sampled)

# # y = [(dNdts_sampled[i_s] - dNdt_complete) / dNdt_complete for i_s in range(N_s)]
# # y = [sum(y[:i]) / i for i in range(1, N_s)]
# # y = [(i**2).sum()**.5 for i in y]

# # plt.plot(y)
# # plt.show()
# # plt.close()

# # y = [(dNdt_sampled - dNdt_complete) / dNdt_complete for i_s in range(N_s)]


# # N_s = 10
# # drhodts_sampled = []
# # for i_s in tqdm(range(N_s)):
# #     kernel = SampledKernel(cfg, rho_complete, R_coag, R_frag)

# #     rho_sampled = solve_react_0d_equation(rho_complete, s, 
# #         rmat=rmat, kmat=kernel.K, dt=dt, niter=N_iter, method="radau"
# #     )
# #     rho_sampled = fix_negative_densities(mg.bin_centers, rho_sampled)

# #     drho_sampled = rho_sampled - rho_complete
# #     drhodt_sampled = drho_sampled / dt
# #     drhodts_sampled.append(drhodt_sampled)

# #     # plt.semilogy(dNdt_sampled)
# #     # plt.show()
# #     # plt.close()


# # y = drhodts_sampled
# # y = [ sum(y[:i]) / i for i in range(1, N_s)]
# # print(y[0])
# # y = [ (i**2).sum()**.5 for i in y]
# # print(y)

# # y = [ (y[i_s]**2).sum()**.5 for i_s in range(N_s) ]
# # plt.plot(y)
# # plt.show()
# # plt.close()




