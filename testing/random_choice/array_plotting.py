import matplotlib.pyplot as plt
from matplotlib import colors


FIGSIZE = [17, 4]


def format_coord(i, j):
    i = int(i)
    j = int(j)
    value = i * j
    text = f"{i=} {j=} {value=}"
    return text


def plot_array(arr, cmap_bounds):
    plt.figure(figsize=FIGSIZE)
    N = len(arr.shape)
    for k in range(N):
        arr_k = arr[k]
        ax = plt.subplot(1, N, k+1)
        ax.set_aspect('equal', adjustable='box')
        ax.format_coord = format_coord
        plt.set_cmap("Reds")
        vmin, vmax = cmap_bounds
        cmap_norm = colors.Normalize(vmin=vmin, vmax=vmax)
        plt.pcolor(arr_k, norm=cmap_norm)
        plt.colorbar()
    plt.show()
    plt.close()
