import mc_coag
from mc_coag.kernel import Kernel

kernel = Kernel()
print(kernel)
K = kernel.entries()
print(K)

kernel = Kernel.from_config()
print(kernel)
