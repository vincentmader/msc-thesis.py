import numpy as np
import matplotlib.pyplot as plt

# Define x-axis.
nx = 1000
x0 = -4
x1 = 4
x  = np.linspace(x0, x1, nx)
dx = x[1] - x[0]

# Define function on y-axis.
f = np.exp(-x**2/2)

# Calculate total integral.
tot = f.sum() * dx
print(tot)

# Define sampling probability.
# p = f / tot * dx
p = np.ones(shape=[nx]) / nx

p = np.abs(x)
# p = p.max() - p
p /= p.sum()
i = np.arange(nx)

# Sample.
nsam = 100
indxs = np.random.choice(i, p=p, size=nsam, replace=True)

# Calculate sampled integral.
tots = (f[indxs] / p[indxs]).sum() * dx / nsam
print(tots)

print(tots / tot)

plt.plot(x, f)
plt.plot(x, p/p.max())
plt.plot(x[indxs], f[indxs], ".")
plt.show()
plt.close()
