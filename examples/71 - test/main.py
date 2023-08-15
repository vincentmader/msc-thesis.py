import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.vX_2023_08_14.base import *
except ModuleNotFoundError as e:
    raise e


cfg = Config()
kernel_1 = Kernel(cfg)
kernel_2 = Kernel(cfg)
kernel_3 = Kernel(cfg)

s1 = KernelSubplot(
    kernel_1, 20, 
    title="pcolor plot 1",
    xlabel="particle radius $a_i$",
    ylabel="particle radius $a_j$",
)
s2 = KernelSubplot(
    kernel_2, 25, 
    title="pcolor plot 2",
    xlabel="particle radius $a_i$",
)
s3 = KernelSubplot(
    kernel_3, 40, 
    title="pcolor plot 3",
    xlabel="particle radius $a_i$",
)
subplots = [
    s1, 
    s2,
    s3,
]

p = GridspecPlot(subplots)
p.render()
