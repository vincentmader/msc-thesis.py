import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from functions.plotting.preset.p1 import main as p1
except ModuleNotFoundError as e:
    raise e

rho_sample = 0.5

cfg = Config(
    # initial_mass_bin=40,
    # mass_resolution=100,
    mass_resolution=50,
    # enable_collision_sampling=False,
    # enable_collision_sampling=False,
    # mass_resolution=100,
)
N_m = cfg.mass_resolution
nr_of_samples = int((N_m**2 + N_m) / 2 * rho_sample)
cfg.nr_of_samples = nr_of_samples


if __name__ == "__main__":
    p1(cfg)
