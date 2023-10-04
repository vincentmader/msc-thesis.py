import numpy as np
from constants import SECONDS_PER_YEAR

# ═════════════════════════════════════════════════════════════════════════════

MILLION, YEAR = 1e6, SECONDS_PER_YEAR
t_max = 100 * MILLION * YEAR
print(f"t_max = {t_max / YEAR:.2e} y")
print(f"      = {t_max:.2e} s")

# ═════════════════════════════════════════════════════════════════════════════

a = 60
b = 70

t_min = 1
t_max = 25

f = (t_max / t_min) ** (1 / 100)

t = np.logspace(t_min, t_max, 100)
print(t)
print()


years = [60, 70]
for year in years:
    duration_in_years = t[year] / SECONDS_PER_YEAR
    duration_in_myears = duration_in_years / 1e6

    print(f"\nt[{year}]")
    print(f"{duration_in_years = :.2E} y")
    print(f"{duration_in_myears = :.2E} My")
