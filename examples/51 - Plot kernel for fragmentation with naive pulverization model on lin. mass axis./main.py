import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from functions.plotting.preset.p1 import main as p1
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=51,
    mass_resolution=50,
    initial_mass_bin=10,

    enable_coagulation=False,
    enable_fragmentation=True,
    enable_cancellation_handling=True,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    fragmentation_variant="naive/pulverization",
    enable_collision_sampling=False,
)

if __name__ == "__main__":
    p1(cfg)
