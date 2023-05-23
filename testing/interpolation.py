import matplotlib.pyplot as plt
import numpy as np

def f(x): return 1 / x

X_MIN, X_MAX, N_X = 1, 2, 100
X = np.logspace(X_MIN, X_MAX, N_X)
Y = f(X)

x = 15
y = np.interp(x, X, Y)

print(y)
print(f(x))

plt.plot(X, Y)
plt.scatter([x], [y])
plt.show()
