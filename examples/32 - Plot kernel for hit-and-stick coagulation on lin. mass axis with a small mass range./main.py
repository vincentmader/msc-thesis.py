import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from functions.plotting.preset.p1 import main as p1
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    # Define mass axis.
    # On a linear grid, if we want to reach a mass grid spacing of exactly one,
    # we have to chose `mass_resolution = mass_max_value - mass_min_value`.
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=11,
    mass_resolution=10,
    # Define processes to include in the simulation.
    enable_coagulation=True,
    enable_fragmentation=False,
    enable_cancellation_handling=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    enable_collision_sampling=False,
)

if __name__ == "__main__":
    p1(cfg)
