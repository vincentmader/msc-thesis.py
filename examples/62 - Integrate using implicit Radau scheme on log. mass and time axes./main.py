import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis, DiscreteRadialAxis, DiscreteTimeAxis
    from config import Config
    from disk import mass_distribution
    from disk.disk import disk_mass_from_distribution
    from kernel import Kernel
    from solver import Solver
    from visualization.evolution import EvolutionPlot, MassConservationPlot
except ModuleNotFoundError as e:
    raise e

# Load configuration from `../../config.toml`.
cfg = Config()

# Define discrete axis for radial distance from star, as well as for mass.
rg = DiscreteRadialAxis(cfg)
mg = DiscreteMassAxis(cfg)
mc = mg.bin_centers
dm = mg.bin_widths

# Define kernel.
kernel = Kernel(cfg)
K = kernel.K

# Define temporal domain & solver.
tg = DiscreteTimeAxis(cfg)

# ─────────────────────────────────────────────────────────────────────────────
# from disk.disk import Disk
# from disk.disk_region import DiskRegion
# from collision import collision_rate
# from kees_kernel import create_coag_kernel
# disk = Disk(cfg, rg, mg)
# disk_region = DiskRegion(cfg, disk)
# Cij = collision_rate(cfg, disk, disk_region)
# mgrain = mg.bin_centers
# K = create_coag_kernel(mgrain, Cij)  # Kees
# ─────────────────────────────────────────────────────────────────────────────


def plot_1(kernel, N, f, m2f, dm2f):
    p = EvolutionPlot(kernel, N, f, m2f, dm2f)
    p.render()


def plot_2(t, Ms):
    p = MassConservationPlot(cfg, t, Ms)
    p.render()


if __name__ == "__main__":

    # Initialize mass distribution.
    n0 = mass_distribution.dirac_delta(cfg)
    # n0 = mass_distribution.mrn_distribution(cfg)

    solver = Solver(cfg)

    # Run the solver.
    N, f, m2f, dm2f = solver.run(n0, K)

    # Prepare abscissa & ordinate for plot of disk mass error.
    t = tg.bin_centers
    Ms = [disk_mass_from_distribution(n, mc, dm) for n in f]

    # Create plots.
    plot_1(kernel, N, f, m2f, dm2f)
    plot_2(t, Ms)
