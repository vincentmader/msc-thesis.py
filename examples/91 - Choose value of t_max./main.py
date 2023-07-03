SECOND = 1
MINUTE = 60 * SECOND
HOUR = 60 * MINUTE
DAY = 24 * HOUR
YEAR = 365.25 * DAY

MILLION = 1e6

t_max = 100 * MILLION * YEAR

print(f"t_max = {t_max / YEAR:.2e} y")
print(f"      = {t_max:.2e} s")
