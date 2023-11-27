import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from functions.plotting.preset import p1
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    initial_mass_bin=10,
    enable_coagulation=False,
    enable_fragmentation=True,
    enable_cancellation_handling=True,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    fragmentation_variant="mrn",
)

if __name__ == "__main__":
    p1(cfg)
