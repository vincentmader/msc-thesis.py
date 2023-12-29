import matplotlib.pyplot as plt
import numpy as np

# Define the matrix dimensions.
N_m = 10
# Create a meshgrid.
i      = np.arange(N_m)
j      = np.arange(N_m)
j      = j
jj, ii = np.meshgrid(i, j)

yy, xx = ii, jj

def main():
    K = heaviside_theta(ii - jj)
    # K *= (ii*2 + jj*2)

    for b in range(3):
        zz = K
        for _ in range(b):
            zz = symmetrize_matrix(zz)
    
        plt.figure(figsize=(8, 8))
        plt.gca().set_aspect('equal', "box")  # = plot aspect ratio
        if N_m < 20:
            plt.xticks(j)
            plt.yticks(i)
        draw_grid()
        
        title = ""
        if b == 0:
            title = r"$K_{kij}=\theta_{H}(i-j)$"
        elif b == 1:
            title = "$X^1:=(K_{kij}+K+{kij}^T)/2$"
        elif b == 2:
            title = "X^2:=$(X^1_{kij}+X^1+{kij}^T)/2$"
        plt.title(title)

        plt.xlabel("bin $j$")
        plt.ylabel("bin $i$")

        plt.imshow(zz, cmap='Blues', origin="lower")
        plt.colorbar()

        plt.show()
        plt.close()

def heaviside_theta(x: np.ndarray):
    out = np.zeros(x.shape)
    out[np.where(x  > 0)] = 1
    out[np.where(x  < 0)] = 0
    out[np.where(x == 0)] = .5
    return out

def symmetrize_matrix(m: np.ndarray):
    return (m + m.T) / 2

def draw_grid():
    for i in range(N_m+1):
        for j in range(N_m+1):
            x = i - 0.5
            y = i - 0.5
            plt.plot([-0.5, N_m], [y, y], 'k--', linewidth=.5)
            plt.plot([x, x], [-0.5, N_m], 'k--', linewidth=.5)

if __name__ == "__main__":
    main()
