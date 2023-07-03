import numpy as np
from numpy import pi as PI

from constants import m_p, k_B, alpha
from disk.dust_particle import particle_radius_from_mass


def dv_turbulence(cfg, disk, disk_region):
    mg = disk.mass_axis
    mc = mg.grid_cell_centers
    radii = particle_radius_from_mass(mc)

    T_mid = disk_region.T_mid
    rho_g = disk_region.rho_g
    t_stop = disk_region.stopping_time(radii)
    c_s = disk_region.c_s
    H_p = disk_region.H_p
    Omega_K = disk_region.Omega_K

    vturb0 = np.sqrt(alpha) * c_s
    lturb0 = np.sqrt(alpha) * H_p
    tturb0 = 1 / Omega_K  # lturb0 / vturb0

    def v1():
        return np.zeros(shape=[mg.N] * 2)

    def v2():
        return dv_tu_voelk(t_stop, t_stop, rho_g, T_mid, vn, tn)

    def v3():
        dv = np.zeros(shape=[mg.N] * 2)
        for i, _ in enumerate(mc):
            for j, _ in enumerate(mc):
                t_i = t_stop[i]
                t_j = t_stop[j]
                dv[i, j] = ormelcuzzi(t_i, t_j, rho_g, T_mid, vturb0, tturb0)
        return dv

    dv = v3()
    return dv


def dv_tu_voelk(ts1, ts2, rhogas, tgas, vturb0, tturb0):
    """
    FITTING FUNCTION TO VOELK TURBULENT DELTAV
          function written by C. Dominik

    Derived from Voelk et al 1980, but using the
    analytical fits by Weidenschilling 1984 and the
    adaptation for small grain by by Mizuno et al.

    ARGUMENTS:
     ts1      Stopping time of particle 1 [s]
     ts2      Stopping time of particle 2 [s]
     rho      Gas density, [g cm^-3]
     tgas     Gas temperature
     vturb0   Turbulent velocity of larges eddie
     tturb0   Turnover time of largest eddie

    RETURNS:
     voelk    The average delta V [cm/s] between the grains

    """
    mu = 2.3
    nH2 = rhogas / mu / m_p
    sig = 5.476e-17 * 4.e0 * PI    # 5.476e-17 = 74d-10**2
    lamH2 = 1.0 / (sig * nH2)
    vtherm = np.sqrt(k_B * tgas / mu / m_p)
    nu = vtherm * lamH2
    Re = vturb0**2 * tturb0 / nu    # Reynolds number

    # Calculate the eddie properties at the smallest scale, where the
    # flow becomes laminar
    vturb_s = vturb0 / np.sqrt(np.sqrt(Re))
    tturb_s = tturb0 / np.sqrt(Re)

    # Sort the stopping times
    t1 = np.array(ts1)
    t2 = np.array(ts2)
    mask = t1 > t2
    t = t1.copy()
    t1[mask] = t2[mask]
    t2[mask] = t[mask]

    # Basic formula
    voelk = vturb0 * 3.0 * t2 / (t1 + t2) * np.sqrt(t2 / tturb0)

    # Adjust if both grains have stopping times smaller
    # than the smallest eddies
    mask = np.logical_and(t1 < tturb_s, t2 < tturb_s)
    voelk[mask] = vturb_s / tturb_s * \
        np.abs(t1 - t2) * np.sqrt(np.log(Re) /
                                  (2.0 * np.sqrt(Re)) * tturb0 / (t1 + t2))
    mask = np.logical_and(mask, Re < 100)
    voelk[mask] = 0.

    # If one large, one small
    mask = np.logical_and(t1 < tturb0, t2 >= tturb0)
    voelk[mask] = vturb0

    # Both large
    mask = np.logical_and(t1 >= tturb0, t2 >= tturb0)
    voelk[mask] = vturb0 * tturb0 * (t1 + t2) / (2.0 * t1 * t2)

    # Limit
    mask = voelk > vturb0
    voelk[mask] = vturb0

    # Check and return
    assert voelk.min() >= 0, 'Error: Turbulent dv negative'
    return voelk


def ormelcuzzi(t1, t2, rho, tgas, vn, tn):
    """
          TURBULENT RELATIVE VELOCITIES ORMEL & CUZZI 2007
                   function written by F. Brauer
                     adapted by C.P. Dullemond

     ARGUMENTS:
      t1       Stopping time of particle 1 [s]
      t2       Stopping time of particle 2 [s]
      rho      Gas density, [g cm^-3]
      tgas     Gas temperature
      vn       Turbulent velocity of largest eddy
      tn       Turnover time of largest eddy (In OC2007 this is called t_L)

     RETURNS:
      ormel    The average delta V [cm/s] between the grains
    """
    mu = 2.3
    mp = 1.6726e-24
    kk = 1.3807e-16
    cs = np.sqrt(kk * tgas / mu / mp)    # Isothermal sound speed
    mmol = mu * mp                     # Mass of H2 molecule
    nmol = rho / mmol                # Nr of H2 molecules / cm^3
    sig = 2e-15                     # Cross section for H2-H2 collisions
    lmfp = 1. / (nmol * sig)             # Mean free path
    vth = np.sqrt(8 / np.pi) * cs     # Thermal velocity of gas particles
    nu = vth * lmfp / 2.0          # Kinematic molecular viscosity
    ln = vn * tn                     # Large eddy scale
    Re = ln * vn / nu                  # Reynolds number
    # Kolmogorov scale (In OC2007 this is called eta)
    ls = Re**(-3. / 4.) * ln
    # Kolmogorov time (In OC2007 this is called t_eta)
    ts = Re**(-0.5) * tn
    vs = ls / ts
    st1 = max(t1, t2) / tn             # Stokes number of largest particle
    tsm = st1 * tn                    # Stopping time of largest particle
    vg = vn * np.sqrt(1.5)   # see ormel paper vg <-> vn
    if (tsm <= ts):
        # Limit of tightly coupled particles (below the Kolmogorov scale)
        # Section 3.4.1 of OC2007, their Eq. 26
        du1 = (t1 / tn)**2 / (t1 / tn + Re**(-0.5))
        du2 = (t2 / tn)**2 / (t2 / tn + Re**(-0.5))
        dv = vg * np.sqrt((t1 - t2) / (t1 + t2) * (du1 - du2))
    else:
        # At least one of the particles is not in the tightly
        # coupled regime.
        if (tsm <= tn):
            # Both particles have St<1
            # Section 3.4.2 of OC2007, their Eq. 28
            # Ratio of St_smallest/St_largest
            eps = min(t1, t2) / max(t1, t2)
            st1 = max(t1, t2) / tn                   # St1 = St_largest
            du1 = 1. / (1. + 1.6) + eps**3 / (1.6 + eps)
            du2 = 3.2 - (1 + eps) + 2 / (1 + eps) * du1
            dv = vg * np.sqrt(du2 * st1)
        else:
            # At least one of the particles has St>=1
            # Section 3.4.3 of OC2007, their Eq. 29
            dv = (tn / (tn + t1) + tn / (tn + t2))**0.5 * vg
    return dv


def ormelcuzzi_arraywise(t1, t2, rho, tgas, vn, tn):
    """
          TURBULENT RELATIVE VELOCITIES ORMEL & CUZZI 2007
                   function written by F. Brauer
                     adapted by C.P. Dullemond

     ARGUMENTS:
      t1       Stopping time of particle 1 [s]
      t2       Stopping time of particle 2 [s]
      rho      Gas density, [g cm^-3]
      tgas     Gas temperature
      vn       Turbulent velocity of largest eddy
      tn       Turnover time of largest eddy (In OC2007 this is called t_L)

     RETURNS:
      ormel    The average delta V [cm/s] between the grains
    """
    mu = 2.3
    mp = 1.6726e-24
    kk = 1.3807e-16
    cs = np.sqrt(kk * tgas / mu / mp)    # Isothermal sound speed
    mmol = mu * mp                     # Mass of H2 molecule
    nmol = rho / mmol                # Nr of H2 molecules / cm^3
    sig = 2e-15                     # Cross section for H2-H2 collisions
    lmfp = 1. / (nmol * sig)             # Mean free path
    vth = np.sqrt(8 / np.pi) * cs     # Thermal velocity of gas particles
    nu = vth * lmfp / 2.0          # Kinematic molecular viscosity
    ln = vn * tn                     # Large eddy scale
    Re = ln * vn / nu                  # Reynolds number
    # Kolmogorov scale (In OC2007 this is called eta)
    ls = Re**(-3. / 4.) * ln
    # Kolmogorov time (In OC2007 this is called t_eta)
    ts = Re**(-0.5) * tn
    vs = ls / ts
    st1 = np.maximum(t1, t2) / tn      # Stokes number of largest particle
    tsm = st1 * tn                    # Stopping time of largest particle
    vg = vn * np.sqrt(1.5)           # see ormel paper vg <-> vn
    # Limit of tightly coupled particles (below the Kolmogorov scale)
    # Section 3.4.1 of OC2007, their Eq. 26
    du1 = (t1 / tn)**2 / (t1 / tn + Re**(-0.5))
    du2 = (t2 / tn)**2 / (t2 / tn + Re**(-0.5))
    dv_tc = vg * np.sqrt((t1 - t2) / (t1 + t2) * (du1 - du2))
    # Both particles have St<1
    # Section 3.4.2 of OC2007, their Eq. 28
    # Ratio of St_smallest/St_largest
    eps = np.minimum(t1, t2) / np.maximum(t1, t2)
    st1 = np.maximum(t1, t2) / tn                     # St1 = St_largest
    du1 = 1. / (1. + 1.6) + eps**3 / (1.6 + eps)
    du2 = 3.2 - (1 + eps) + 2 / (1 + eps) * du1
    dv_sm = vg * np.sqrt(du2 * st1)
    # At least one of the particles has St>=1
    # Section 3.4.3 of OC2007, their Eq. 29
    dv_bg = (tn / (tn + t1) + tn / (tn + t2))**0.5 * vg
    # Combine
    dv = dv_tc.copy()
    dv[tsm > ts] = dv_sm[tsm > ts]
    dv[tsm > tn] = dv_bg[tsm > tn]
    return dv


def get_ormelcuzzi_for_disk(mstar_msun=1, r_au=1, sigma_g=100, hpr=0.025, alpha=1e-2, Stmin=1e-8, Stmax=1e2, nSt=20):
    GG = 6.67408e-08    # Gravitational constant  [cm^3/g/s^2]
    mp = 1.6726e-24     # Mass of proton          [g]
    kk = 1.3807e-16     # Bolzmann's constant     [erg/K]
    MS = 1.98892e33     # Solar mass              [g]
    au = 1.49598e13     # Astronomical Unit       [cm]
    mstar = mstar_msun * MS
    r = r_au * au
    omk = np.sqrt(GG * mstar / r**3)
    hp = hpr * r
    cs = hp * omk
    tgas = cs**2 * 2.3 * mp / kk
    rho = sigma_g / (np.sqrt(2 * np.pi) * hp)
    tn = 1 / omk
    vn = np.sqrt(alpha) * cs
    St = Stmin * (Stmax / Stmin)**np.linspace(0, 1, nSt)
    zeros = np.zeros(nSt)
    St1 = St[:, None] + zeros[None, :]
    St2 = St[None, :] + zeros[:, None]
    t1 = St1 / omk
    t2 = St2 / omk
    dv = ormelcuzzi_arraywise(t1, t2, rho, tgas, vn, tn)
    return dv, St1, St2


def plot_ormelcuzzi_for_disk(mstar_msun=1, r_au=1, sigma_g=100, hpr=0.025, alpha=1e-2, Stmin=1e-8, Stmax=1e2, nSt=100):
    import matplotlib.pyplot as plt
    import numpy as np
    dv, St1, St2 = get_ormelcuzzi_for_disk(
        mstar_msun=mstar_msun, r_au=r_au, sigma_g=sigma_g, hpr=hpr, alpha=alpha, Stmin=Stmin, Stmax=Stmax, nSt=nSt)
    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    plt.contour(np.log10(St1), np.log10(St2), np.log10(dv), levels=20)
    plt.xlabel(r'$^{10}\log(\mathrm{St}_1)$')
    plt.ylabel(r'$^{10}\log(\mathrm{St}_2)$')
