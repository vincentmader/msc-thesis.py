#!../../.venv/bin/python3
import os, sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteMassAxis
    from axis import AxisLabelVariant
    from config import Config
    from kernel import Kernel
    from visualization.preset import p1
except ModuleNotFoundError as e:
    raise e

cfg = Config(
    enable_collision_sampling=True,
    enable_coagulation=True,
    enable_fragmentation=True,
    initial_mass_bin=0,
    # enable_physical_collisions=False,
)

mg, kernel = DiscreteMassAxis(cfg), Kernel(cfg)
scale = mg.scale
axis_label_variant = AxisLabelVariant.Radius if scale == "log" else AxisLabelVariant.Bin
z_limits = (1e-20, 1e+10) if scale == "log" else (-1, 1)
R = kernel.R_coag + kernel.R_frag

def _integrate():
    global t, f, N, m2f, dm2fdt, M
    t, f, N, m2f, dm2fdt, M = p1.integrate(cfg, kernel)

def plot_total_kernel():
    p1.plot_kernel(cfg, mg, kernel, scale, z_limits, axis_label_variant)

def plot_kernel_gain_loss():
    # Plot K_gain & K_loss  with log. colorscale.
    p1.plot_kernel_gain_loss(cfg, mg, kernel, scale, z_limits, axis_label_variant=axis_label_variant)

def plot_kernel_error():
    p1.plot_kernel_error(cfg, mg, kernel, R, scale, z_limits, axis_label_variant=axis_label_variant)

def plot_evolution():
    p1.plot_evolution(cfg, mg, kernel, scale, t, N, f, m2f, dm2fdt)
    # plot_surface(cfg, mg, kernel, scale, t, f, N, m2f, dm2f)

def plot_error():
    p1.plot_error(cfg, mg, kernel, t, M)

def toggle_sampling():
    cfg.enable_collision_sampling = not cfg.enable_collision_sampling

# plot_total_kernel,
plot_kernel_gain_loss()
plot_kernel_error()
_integrate()
plot_evolution()
plot_error()
