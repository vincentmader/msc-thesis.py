from config import Config
from kernel import Kernel
from visualization.v1.mass_conservation import KernelMassConservationPlot

cfg = Config()
kernel_1 = Kernel(cfg)
K = kernel_1.K

mg = kernel_1.mg

def main():
    # Plot `\sum_{ij} m_k \Delta m_k K_kij`.
    p = KernelMassConservationPlot(cfg, mg, K)
    p.show()
