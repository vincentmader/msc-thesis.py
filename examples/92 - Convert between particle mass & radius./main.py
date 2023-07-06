import os
import sys
try:
    sys.path.append(os.path.join("..", "..", "src"))
    from config import Config
    from dust import particle_mass_from_radius
except ModuleNotFoundError as e:
    raise e

centimeter = 0.01

if __name__ == "__main__":
    cfg = Config()

    r_min = 1e-4 * centimeter
    r_max = 1e+3 * centimeter

    print(f"r_min = {r_min/centimeter:.2e} cm")
    print(f"r_max = {r_max/centimeter:.2e} cm\n")

    print(f"r_min = {r_min:.2e} m")
    print(f"r_max = {r_max:.2e} m\n")

    m_min = particle_mass_from_radius(r_min)
    m_max = particle_mass_from_radius(r_max)
    print(f"m_min = {m_min:.2e} kg")
    print(f"m_max = {m_max:.2e} kg")

    # m_min = cfg.mass_min_value
    # m_max = cfg.mass_max_value
    # r_min = particle_radius_from_mass(m_min)
    # r_max = particle_radius_from_mass(m_max)
    # print(f"r_min = {r_min:.2e}")
    # print(f"r_max = {r_max:.2e}")
