import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from visualization.preset import p1
except ModuleNotFoundError as e:
    raise e

cfg = Config()


if __name__ == "__main__":
    p1(cfg)
