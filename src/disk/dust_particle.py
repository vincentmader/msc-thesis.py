from numpy import pi as PI

from constants import rho_s


def particle_radius_from_mass(m):  # TODO Redefine m_min & m_max
    V_m = m / rho_s
    a_m = (V_m * 3 / 4 / PI)**(1 / 3)
    return a_m


def particle_mass_from_radius(a):
    V_a = 4 / 3 * PI * a**3
    m_a = rho_s * V_a
    return m_a
