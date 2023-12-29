import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

X_MIN, X_MAX      = (-5, +5)
NR_OF_X_VALUES    = 1000
NR_OF_X_SAMPLES   = 100
NR_OF_SAMPLE_RUNS = 1000

def probability_from_homogenous_distribution(x):
    y = np.ones(shape=len(x))
    return y / y.sum()
def probability_from_function(x, f):
    y = f(x)
    return y / y.sum()
def probability_from_triangle_good(x):
    y = np.abs(x)
    y = y.max() - y
    return y / y.sum()
def probability_from_triangle_bad(x):
    y = np.abs(x)
    return y / y.sum()

FUNCTIONS = {
    # "exponential"   : lambda x: np.exp(x),
    "gaussian"      : lambda x: np.exp(-x**2 / 2),
}

SAMPLING_PROBABILITIES = {
    "homogenous"    : probability_from_homogenous_distribution,
    "function"      : probability_from_function,
    "triangle_good" : probability_from_triangle_good,
    "triangle_bad"  : probability_from_triangle_bad,
}

x  = np.linspace(X_MIN, X_MAX, NR_OF_X_VALUES)
dx = x[1] - x[0]  # Note: This only works because we're on a linear grid.

mu_vs_title = {}
sigma_vs_title = {}

for function_name, function in tqdm(FUNCTIONS.items()):
    y = function(x)
    A_total = np.sum(y * dx)

    for probability_name, probability in SAMPLING_PROBABILITIES.items():
        for replace in [True]:  # False
            title = f"function={function_name}, probability={probability_name}, {replace=}"
            print(title)

            p = probability(x) if probability_name != "function" else probability(x, function)

            A_sampled_vs_sample_run = []
            for _ in range(NR_OF_SAMPLE_RUNS):
                indices_all     = np.arange(NR_OF_X_VALUES)
                indices_sampled = np.random.choice(indices_all, p=p, size=NR_OF_X_SAMPLES, replace=replace)

                x_sampled       = x[indices_sampled]
                y_sampled       = y[indices_sampled]
                p_sampled       = p[indices_sampled]

                A_sampled   = (y_sampled * dx / p_sampled).sum() / NR_OF_X_SAMPLES 
                A_sampled_vs_sample_run.append(A_sampled)      # ^ TODO Why divide by `NR_OF_SAMPLES` here?
            A_sampled_vs_sample_run = np.array(A_sampled_vs_sample_run)

            # ═════════════════════════════════════════════════════════════════════════════

            x_plot = x
            y_plot = y
            plt.plot(x_plot, y_plot, label="$f(x)$")

            x_plot = x
            y_plot = p / p.max()
            plt.plot(x_plot, y_plot, label="$p(x) / p_{max}$")
            
            x_plot = x_sampled
            y_plot = y_sampled
            plt.plot(x_plot, y_plot, '.', label="samples")

            plt.title(title)
            plt.xlabel("$x$")
            plt.ylabel("$y$")
            plt.legend(loc="upper right")
            plt.show()
            plt.close()

            # ═════════════════════════════════════════════════════════════════════════════
            
            x_plot = np.arange(NR_OF_SAMPLE_RUNS)
            y_plot = np.array([A_total] * NR_OF_SAMPLE_RUNS)
            plt.plot(x_plot, y_plot, label="$A_{total}$    $=$" + f"{A_total}")

            A_mean = np.mean(A_sampled_vs_sample_run)
            sigma  = ( np.sum( (A_sampled_vs_sample_run - A_mean)**2 ) / NR_OF_SAMPLE_RUNS )**.5

            mu_vs_title[title] = A_mean
            sigma_vs_title[title] = sigma
            
            x_plot = np.arange(NR_OF_SAMPLE_RUNS)
            y_plot = A_sampled_vs_sample_run # - A_total
            plt.plot(x_plot, y_plot, '.', label="$A_{sampled}=$" + f"{A_mean:.2g} +- {sigma:.2g}")

            plt.title(title)
            plt.xlabel("index of sampling run")
            plt.ylabel("result of sampling run: integral $A_{sampled}$")
            plt.legend(loc="upper right")
            plt.show()
            plt.close()

    # ═════════════════════════════════════════════════════════════════════════════

    titles  = list(sigma_vs_title.keys())
    sigmas  = [sigma_vs_title[t] for t in titles]
    mus     = [   mu_vs_title[t] for t in titles]
    x_plot  = np.arange(len(titles))

    plt.figure(figsize=(16, 5))

    plt.subplot(1, 3, 1)
    plt.title(r"$\mu(A_{sampled,i})$ (mean)")
    plt.scatter(mus, x_plot)

    plt.yticks(x_plot, titles)

    plt.subplot(1, 3, 2)
    plt.title(r"$\sigma(A_{sampled,i})$ (std. dev.)")
    plt.scatter(sigmas, x_plot)

    plt.subplot(1, 3, 3)
    plt.title(r"$\frac{\mu(A_{sampled,i})-A_{total}}{A_{total}}$ (rel. dev.)")
    plt.scatter((mus - A_total) / A_total, x_plot)

    plt.tight_layout()

    plt.show()
    plt.close()
