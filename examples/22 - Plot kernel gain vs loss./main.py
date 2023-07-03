import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.kernel.v2.kernel_gain_vs_loss import KernelGainVsLossPlot
except ModuleNotFoundError as e:
    raise e


cfg = Config(
    enable_fragmentation=True,
    enable_coagulation=True,
)
kernel = Kernel(cfg)

p = KernelGainVsLossPlot(kernel, enable_slider=True)
p.show()
