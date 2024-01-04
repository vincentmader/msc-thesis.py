from dataclasses import field
import sys

import numpy as np
from tqdm import tqdm

from config import Config
from config import PATH_TO_COAG, SOLVERS
from functions.disk import mass_distribution
from functions.dust.collision import collision_outcome_probabilities, collision_rate
from functions.dust.relative_velocity import relative_velocity
from functions.solver import fix_negative_densities
from models.axis import DiscreteMassAxis, DiscreteRadialAxis, DiscreteTimeAxis
from models.disk import Disk, DiskRegion
from models.disk import Disk, DiskRegion
from models.kernel import Kernel, SampledKernel

sys.path.append(PATH_TO_COAG)
from coag.react_0d import solve_react_0d_equation


N_subst = 1  # Nr of time substeps between storage of result
N_iter  = 4  # Nr of iterations for implicit time step


class SolverV2():

    cfg:                    Config
    disk:                   Disk
    region:                 DiskRegion
    rg:                     DiscreteRadialAxis
    mg:                     DiscreteMassAxis
    tg:                     DiscreteTimeAxis
    R_coag:                 np.ndarray
    R_frag:                 np.ndarray
    W_ij:                   np.ndarray

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

        self.cfg            = cfg
        self.disk           = Disk(cfg)
        self.region         = DiskRegion(cfg, self.disk)
        self.rg             = DiscreteRadialAxis(cfg)
        self.mg             = DiscreteMassAxis(cfg)
        self.tg             = DiscreteTimeAxis(cfg)

        # Initialize result matrices with zeroes.
        self.n_k_vs_t       = np.zeros(shape=(self.tg.N, self.mg.N))
        self.N_k_vs_t       = np.zeros(shape=(self.tg.N, self.mg.N))
        self.M_k_vs_t       = np.zeros(shape=(self.tg.N, self.mg.N))
        self.dMdt_vs_t      = np.zeros(shape=(self.tg.N, self.mg.N))
        self.P_ij_vs_t      = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))
        self.S_ij_vs_t      = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))
        self.dK_ij_vs_t     = np.zeros(shape=(self.tg.N, self.mg.N, self.mg.N))

        # Define PPD, & radial position of interest in it.
        disk        = Disk(cfg, self.rg, self.mg) 
        disk_region = DiskRegion(cfg, disk)

        # Define relative dust particle velocities, & collision rates.
        dv          = relative_velocity(cfg, disk, disk_region)
        R_coll      = collision_rate(cfg, disk, disk_region)

        # Define rate of coag. & frag. events.
        P_coag, P_frag = collision_outcome_probabilities(cfg, dv)
        self.R_coag = R_coll * P_coag
        self.R_frag = R_coll * P_frag

        # Define initial state, initial kernel, & sampling weights.
        self.N_k_vs_t[0, :] = mass_distribution.dirac_delta(cfg) * self.mg.bin_widths
        self.kernels        = [Kernel(cfg, self.R_coag, self.R_frag)]
        self.W_ij           = weights(self.mg, self.kernels[0])

    def run(self):

        # Read information about discrete axes from `self`.
        N_m, mc, dm = self.mg.N, self.mg.bin_centers, self.mg.bin_widths
        N_t, tc     = self.tg.N, self.tg.bin_centers

        # Integrate the Smoluchowski equation.
        s = np.zeros((N_m))
        rmat = np.zeros((N_m, N_m))
        for i_t in tqdm(range(1, N_t)):
            self.iteration_idx = i_t
            self.step(rmat, s)

        # Translate units and save results to `self`.
        n    = self.N_k_vs_t / dm
        M    = self.N_k_vs_t * mc
        dMdt = np.array([
            (M[i] - M[i-1]) / (tc[i] - tc[i-1]) for i in range(1, self.tg.N)
        ])
        self.n_k_vs_t       = n
        self.M_k_vs_t       = M
        self.dMdt_vs_t[1:]  = dMdt

    def step(self, rmat, s):

        # Read information about discretized time axis from `self`.
        i_t = self.iteration_idx
        tc = self.tg.bin_centers
        dt = (tc[i_t] - tc[i_t - 1]) / N_subst
        
        # Read current particle density from `self`.
        N_dust = self.N_k_vs_t[i_t - 1]

        # Read kernel from `self` (sample collisions, if enabled).
        if self.cfg.enable_collision_sampling:
            self.sample()
        K = self.kernels[-1].K

        # Loop over substeps.
        for _ in range(N_subst):
            assert self.cfg.solver_variant in SOLVERS, f"{self.cfg.solver_variant} = "

            # Apply explicit Euler scheme.
            if self.cfg.solver_variant == "explicit_euler":
                N_i = N_dust[None, :, None]
                N_j = N_dust[None, None, :]
                dNdt = (K[:, :, :] * N_i * N_j).sum(axis=2).sum(axis=1)
                N_dust += dNdt * dt

            # Apply implicit Euler scheme.
            if self.cfg.solver_variant == "implicit_euler":
                N_dust = solve_react_0d_equation(N_dust, s, 
                    rmat=rmat, kmat=K, dt=dt, niter=N_iter, method="backwardeuler"
                )

            # Apply implicit Radau scheme.
            if self.cfg.solver_variant == "implicit_radau":
                N_dust = solve_react_0d_equation(N_dust, s, 
                    rmat=rmat, kmat=K, dt=dt, niter=N_iter, method="radau"
                )

        # Fix sub-zero densities, and save results to state-vector in `self`.
        N_dust = fix_negative_densities(self.mg.bin_centers, N_dust)
        self.N_k_vs_t[i_t, :] = N_dust.copy()

    def sample(self):

        i_t = self.iteration_idx
        sampled = SampledKernel(
            self.cfg, 
            self.N_k_vs_t[i_t-1], 
            self.R_coag, 
            self.R_frag, 
            W_ij=self.W_ij
        )
        self.kernels.append(sampled)
        self.P_ij_vs_t[i_t] = sampled.P_ij
        self.S_ij_vs_t[i_t] = sampled.N_ij  # TODO Rename.


def weights(mg: DiscreteMassAxis, kernel: Kernel):
    return np.sum([
        mg.bin_centers[k] * np.abs(kernel.K[k]) for k in range(mg.N)
    ])


