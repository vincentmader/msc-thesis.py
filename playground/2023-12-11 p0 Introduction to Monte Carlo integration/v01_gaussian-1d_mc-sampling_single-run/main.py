import matplotlib.pyplot as plt
import numpy as np

# Define x-axis
# ═════════════════════════════════════════════════════════════════════════════
N_x   = 1000                                                             # NOTE
x_min = -5                                                               # NOTE
x_max = +5                                                               # NOTE
x     = np.linspace(x_min, x_max, N_x)
dx    = x[1] - x[0]  # Note: This only works because we're on a linear grid.

# Define y-axis (Gaussian function)
# ═════════════════════════════════════════════════════════════════════════════
f_G   = lambda x: np.exp(-x**2 / 2)
y     = f_G(x)

# Calculate total integral (without sampling) by summing f over entire x-axis.
# ═════════════════════════════════════════════════════════════════════════════
A_tot = np.sum(y * dx)

# Define sampling probabilities.
# ═════════════════════════════════════════════════════════════════════════════
def sampling_probability(variant: str):
    if variant == "function":
        p = y / A_tot * dx
    elif variant == "homogenous":
        p = np.ones(shape=[N_x]) / N_x
    elif variant == "triangle-good":
        p = np.abs(x)
        p = p.max() - p
        p /= p.sum()
    elif variant == "triangle-bad":
        p = np.abs(x)
        p /= p.sum()
    else:
        raise Exception()
    return p

# Perform sampling of x-values.
# ═════════════════════════════════════════════════════════════════════════════
variant = "homogenous"                                                   # NOTE 
replace = True                                                           # NOTE
p = sampling_probability(variant)
i_x   = np.arange(N_x)
N_s   = 100                                                              # NOTE
i_s   = np.random.choice(i_x, p=p, size=N_s, replace=replace)
x_s   = x[i_s]
y_s   = y[i_s]

# Calculate integral from sampled values.
# ═════════════════════════════════════════════════════════════════════════════
A_s   = y[i_s] * dx / N_s
if replace:
    A_s /= p[i_s]
A_s   = A_s.sum()

# Calculate relative deviation/error.
# ═════════════════════════════════════════════════════════════════════════════
err   = (A_s - A_tot) / A_tot

# Print results.
# ═════════════════════════════════════════════════════════════════════════════
print(f"\nA_tot = {A_tot}")
print(f"A_s   = {A_s}")
print(f"\nerr   = (A_s - A_tot) / A_tot\n      = {err}")

# Plot.
# ═════════════════════════════════════════════════════════════════════════════
plt.plot(x, y, label="$f(x)$")
plt.plot(x, p/p.max(), label="$p(x)$")
plt.plot(x_s, y_s, '.', label="samples")

plt.legend()
plt.show()
plt.close
