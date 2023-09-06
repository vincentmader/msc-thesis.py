import os
import sys

import numpy as np
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from kernel import Kernel
except ModuleNotFoundError as e:
    raise e

cfg = Config()
kernel = Kernel(cfg)

# K* = K.abs().normalize()

# P = prob_from_kernel(K)


# (i, j) = random(P)

# ks = affected_ks(i, j)
# for k in ks:
#     K[k, i, j] += R ?





# effect only noticable for large kernel sizes
# 50 ?
# 100 ?
# 200 ?
# how much?





# def prob_dist(n, K):
#     P = [[0]]
#     for i in range(N):
#         for j in range(N):
#             P[i, j] = sum_k K[k,i,j] * n_i * n_j
#     P = P.normalize()
# 
# def pairs(K, p):
#     # K = K.flatten()
#     indices = indices(K)
#     pairs = random.choice(indices, p=P)
#     return pairs
#
# def random_sample_sparse_K():
#     p = prob_dist(n, K)
#     pairs = pairs(K, p)
#     K = [[0]]
#     for (i, j) in pairs:
#         for k in ks(i, j):
#             K[k, i, j] += R

# def integrate(K):
#     for t in ts:
#         K = random_sample_sparse_K()
#         update_mass_distribution(K)


#    P ~ K * n_i * n_j
#    P ~ norm(K) = K / sum(K)

# 1. P ~ K
# 2. P ~ (K_gain + abs(K_loss))

# P = ?
# R = ?

# relevant_collisions = random.choice(P)

# K = 0
# for (i, j) in relevant_collisions:
#     ks = ks(i, j)
#     for k in ks:
#         K[k, i, j] += R

