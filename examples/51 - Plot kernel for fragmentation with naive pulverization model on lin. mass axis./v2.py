import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import KernelAxis
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot
    from visualization.kernel.kernel import KernelSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    mass_axis_scale="lin",
    mass_min_value=1,
    mass_max_value=50,
    enable_coagulation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
    fragmentation_variant="naive/pulverization",
)
kernel = Kernel(cfg)
mg = kernel.mg


def main():
    GridspecPlot([
        KernelSubplot(
            mg, kernel.K_gain, 
            title="kernel gain contribution $G_{kij}$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            axis=KernelAxis.Bin,
        ),
        KernelSubplot(
            mg, -kernel.K_loss,
            title="kernel loss contribution $L_{kij}$",
            scales=("lin", "lin", "lin"),
            symmetrized=True,
            axis=KernelAxis.Bin,
        ),
    ], add_slider=True).render()
