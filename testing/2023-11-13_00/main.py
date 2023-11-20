import numpy as np

N = 50
M = 1000

def v1():
    arr = []
    for i in range(N):
        for j in range(N):
            arr.append((i, j))
    return arr

def v2():
    return [(i, j) for i in range(N) for j in range(N)]

def v3():
    return np.array(np.meshgrid(np.arange(N), np.arange(N))).T.reshape(-1, 2)


assert v1() == v2() == v3()
print("Success!")


import time

for f in [v1, v2, v3]:
    t = 0
    for _ in range(M):
        t_i = time.time()
        arr = f()
        t_f = time.time()
        t += t_f - t_i
    print(t / M)
