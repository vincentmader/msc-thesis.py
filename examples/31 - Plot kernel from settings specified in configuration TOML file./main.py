import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from functions.plotting.preset.p1 import main as p1
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    # initial_mass_bin=40,
)


if __name__ == "__main__":
    p1(cfg)
