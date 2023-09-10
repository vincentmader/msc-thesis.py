import os
import sys

import matplotlib.pyplot as plt
import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from axis import DiscreteTimeAxis
    from config import Config
    from disk import mass_distribution
    from dust import particle_radius_from_mass
    from kernel import Kernel
    from solver import Solver
    from visualization.evolution.v1.slider_plot_2 import InteractiveSliderLinePlot
    from visualization.kernel.v3_2023_08_14.pcolor_matrix_subplot import PcolorMatrixSubplot
    from visualization.kernel.v3_2023_08_14.gridspec_plot import GridspecPlot
except ModuleNotFoundError as e:
    raise e

# First, define the weights matrix: `W_ij = \sum_k m_k abs(K_kij)`.
def _weights(kernel: Kernel) -> np.ndarray:
    r"""Return matrix $W_{ij} = \sum_k m_k |K_{kij}|$."""
    K = kernel.K
    mg = kernel.mg
    mc = mg.grid_cell_centers
    N_m = mg.N

    W_ij = np.zeros(shape=(N_m, N_m))
    for k in range(N_m):
        W_ij += mc[k] * np.abs(K[k])
    return W_ij


# Now, define the probability via normalization: `\sum_{ij} W_ij = 1`
def _probability(weights: np.ndarray) -> np.ndarray:
    r"""Return probability matrix $P_{ij} = W_{ij} / \sum W_{ij}$."""
    return weights / weights.sum()


def _index_matrix(
    shape_2d: tuple[int, int],
) -> np.ndarray:
    """Return 'index matrix', i.e.
        | 0,0 | 1,0 |  ->   | 0 | 1 |
        | 0,1 | 1,1 |  ->   | 2 | 3 |
    """
    N_i, N_j = shape_2d
    I = []
    for i in range(N_i):
        I_i = []
        for j in range(N_j):
            I_ij = i * N_i + j
            I_i.append(I_ij)
        I.append(I_i)
    return np.array(I)


def _sample_indices(
    P: np.ndarray,
    shape_2d: tuple[int, int],
) -> list[int]:
    """Return list of indices corresponding to sampled particle pairs."""
    I = _index_matrix(shape_2d)

    shape_1d = (shape_2d[0] * shape_2d[1])
    P = P.reshape(shape_1d)
    I = I.reshape(shape_1d)

    samples = np.random.choice(I, p=P, size=N_sample)
    return list(samples)


def _ijs_from_indices(
    indices: list[int],
    shape_2d: tuple[int, int],
) -> list[tuple[int, int]]:
    """Return inversely transformed matrix from `self._index_matrix()`.
        [0, 1, 2, 3] -> [(0,0), (1,0), (0,1), (1,1)}]
    """
    out = []
    for idx in indices:
        i = idx // shape_2d[0]
        j = idx % shape_2d[0]
        out.append((i, j))
    return out


def _run_integrator(kernel, K):
    mg = kernel.mg
    solver = Solver(cfg)
    n0 = mass_distribution.dirac_delta(cfg)
    K = K # TODO: Do `K = kernel.K_sample` or something like that.
          # TODO: Write `kernel.K_sample` method.
          #       Maybe in new class `SampledKernel` ?
          #       -> Needs to inherit kernel construction methods.
          #       -> Need to define `K_kij` determination methods (?)
    N, f, m2f = solver.run(mg, n0, K)
    return N, f, m2f


def plot_1(m, m2f, dm2f):
    plot = InteractiveSliderLinePlot(
        cfg,
        m, m2f, dm2f,
        title="Temporal Evolution of Particle Mass Distribution Function",
        ylabel_1="$n_i\cdot m_i\cdot\Delta m_i$    [kg m $^{-3}$]",
        xlabel_1="mass $m_i$ [kg]",
        ylabel_2=r"$\frac{\Delta n_i\cdot m_i\cdot\Delta m_i}{\Delta t}$   [kg m$^{-3}$ s$^{-1}$]",
        xlabel_2="mass $m_i$ [kg]",
        xlims_1=(m[0], m[-1]),
        xlims_2=(m[0], m[-1]),
        ylims_1=[1e-15, 1e-8],
        ylims_2=[1e-40, 1e-15],
    )
    plot.draw()
    plt.show()
    plt.close()


if __name__ == "__main__":

    N_sample = 100

    cfg = Config()
    kernel = Kernel(cfg)

    W = _weights(kernel)
    P = _probability(W)

    # TODO Give `P` to integrator.
    # TODO In each time step, do `P -> P * N_i * N_j` or sth. like that.
    # TODO Call `_sample_indices` and `_ijs_from_indices`.
    # TODO Build kernel from only those `[(i,j), ...]`.
    # TODO Integrate with the new kernel.
    # TODO Compare results to earlier, as well as to analytical solutions.

    shape_2d = P.shape
    indices = _sample_indices(P, shape_2d)
    ijs = _ijs_from_indices(indices, shape_2d)

    kernel = Kernel(cfg, ijs=ijs)
    K = kernel.K
    mc = kernel.mg.grid_cell_centers
    rho_s = cfg.dust_particle_density
    ac = particle_radius_from_mass(mc, rho_s)

    s1 = PcolorMatrixSubplot(
        ac, ac, kernel.K_gain, 
        title="kernel gain contribution $G_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        ylabel="particle radius $a_i$ [m]",
        symmetrize=True,
    )
    s2 = PcolorMatrixSubplot(
        ac, ac, -kernel.K_loss,
        title="kernel loss contribution $L_{kij}$",
        xlabel="particle radius $a_j$ [m]",
        symmetrize=True,
    )
    subplots = [s1, s2]
    p = GridspecPlot(subplots, add_slider=True)
    p.render()

    N, f, m2f = _run_integrator(kernel, K)

    tg = DiscreteTimeAxis(cfg)
    # Calculate temporal derivative of mass distribution.
    dm2f = m2f[1:] - m2f[:-1]
    dm2f = list(dm2f)
    dm2f.append(dm2f[-1])  # TODO Fix array shapes in a better way than this.
    dm2f = np.array(dm2f)
    dm2f = [dm2f[i] / tg.grid_cell_widths[i] for i, _ in enumerate(dm2f)]

    plot_1(mc, m2f, dm2f)
