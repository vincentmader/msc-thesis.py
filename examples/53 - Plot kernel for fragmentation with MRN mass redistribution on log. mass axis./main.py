import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
    from visualization.base import GridspecPlot, PcolorMatrixSubplot
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    # fragmentation_variant="naive/pulverization",
    enable_coagulation=False,
    enable_physical_collisions=False,
    relative_velocity_components=[],
)
kernel_1 = Kernel(cfg)

mg = kernel_1.mg
mc = mg.grid_cell_centers
ac = mg.particle_radii


if __name__ == "__main__":
    GridspecPlot([
        PcolorMatrixSubplot(
            ac, ac, kernel_1.K_gain, 
            title="kernel gain contribution $G_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            ylabel="particle radius $a_i$ [m]",
            symmetrized=True,
        ),
        PcolorMatrixSubplot(
            ac, ac, -kernel_1.K_loss,
            title="kernel loss contribution $L_{kij}$",
            xlabel="particle radius $a_j$ [m]",
            symmetrized=True,
        ),
    ], add_slider=True).render()
