import numpy as np


def create_coag_kernel(mgrain, Cij, fixcancel=True):
    """
    Given a mass grid mgrain and a collision rate matrix Cij, this
    computes the coagulation kernel K_kij for perfect sticking. 

    The percect-sticking coagulation equation (Smoluchowski equation) is

     (d/dt) n(m)  = (1/2) int C(m',m'') delta(m'+m''-m) n(m') n(m'') dm'dm''
                    - int C(m,m') n(m) n(m') dm'

                  = (1/2) int C(m',m'') delta(m'+m''-m) n(m') n(m'') dm'dm''
                        - int C(m,m')   delta(m''-m)    n(m') n(m'') dm'dm''

                  = int K(m,m',m'') n(m') n(m'') dm'dm''

    where C(m',m'') is the collision rate between particles of masses m'
    and m'', and the last line defines the coagulation kernel K(m,m',m'').
    Now we make this discrete, and we replace n(m_k) with N_k / Delta m_k
    where Delta m_k is the mass bin width (see dmgrain array). Now multiply
    both sides by dm:

     (d/dt) n(m)dm = int K(m,m',m'') n(m') n(m'') dm dm'dm''

     (d/dt) N(m)   = int K(m,m',m'') dm N(m') N(m'')

    which, if we define Kkij = K(m_k,m_i,m_j) dm_k, equals

     (d/dt) N_k    = Sum_ij K_kij N_i N_j

    Now Kkij = K(m_k,m_i,m_j) dm_k looks a bit like an ugly definition,
    because of the dm_k on the rhs. But note that the K(m_k,m_i,m_j)
    still has the Dirac delta function in it. We can rewrite the two
    Dirac delta functions times the dm as:

     delta(m_i+m_j-m_k) dm_k = Dkij
     delta(m_j-m_k) dm_k = delta_{j,k}       (delta_{.,.} is the Kronecker delta)

    Now let us, for convenience, link this back to discrete versions
    of the coagulation kernel C(m',m''). Define

     Cij  = C(m_i,m_j) and Ckij = C(m_i,m_j) * Dkij

    Note that Sum_k Ckij = Cij. The Smoluchowski equation in
    discrete form is now

     (d/dt) N_k  = (1/2) Sum_ij Ckij N_i N_j - Sum_i Cik N_i N_k

    Comparing this to the generic equation

     (d/dt) N_k    = Sum_ij K_kij N_i N_j

    we see that

     Kkij = 0.5*Ckij - delta_{j,k}Cik

    This tensor is asymmetric in ij. We can make it symmetric without
    loss of generality:

     Kkij = 0.5*Ckij - 0.5 * ( delta_{j,k}Cik + delta_{i,k}Cjk )

    It is this tensor (=3D array) that is returned.

    ARGUMENTS:

      mgrain       Array of mass grid points

      Cij          2D Array of collision rates Cij = C(m_i,m_j)

      fixcancel    If true, then we fix the problematic near-cancellation
                   occurring if a tiny grain hits a huge grain. The 
                   (1-eps) factor from the gain term will be offset by
                   the 1 factor from the loss term. If eps is too small,
                   then (1-eps)-1 may not become eps, as it should. So
                   if, for these pairs of masses, subtract 1 from both
                   the gain and loss term, the problem is solved. That
                   is done if fixcancel=True. To see the difference,
                   make a mass grid spanning 1e20, and plot:
                     plt.plot(Kkij[-10,:,-10]/mgrain)
                   for the case with fixcancel=False and fixcancel=True.
                   You will see that for fixcancel=False the left side
                   of the plot goes berserk, while for fixcancel=True
                   it stays constant.                   

    RETURNS:

      Kkij         The perfect sticking coagulation kernel
    """
    nm = len(mgrain)
    Dkij, icancel = create_coag_delta(mgrain, fixcancel=fixcancel)
    Kkij = 0.5 * Cij[None, :, :] * Dkij[:, :, :]
    for k in range(nm):
        keep = 1 - icancel[:, k]  # Only [:,k] is correct, not [k,:]!
        Kkij[k, :, k] -= 0.5 * Cij[:, k] * keep
        Kkij[k, k, :] -= 0.5 * Cij[k, :] * keep
    return Kkij


def create_ilow_eps_icancel(mgrain):
    """
    Compute where m_1+m_2 ends up in the mass grid. This is required for the
    perfect sticking coagulation kernel. The ilow is the index of the mass
    point just below m_1+m_2, eps is the linear interpolation coefficient,
    and icancel is a mask saying where ilow(i,j)=i. The icancel is meant to
    handle cases with extreme mass ratios.
    """
    nm = len(mgrain)
    m1, m2 = np.meshgrid(mgrain, mgrain, indexing='ij')
    m = m1 + m2
    msmall = np.zeros((nm, nm))
    msmall[m1 < m2] = m1[m1 < m2]
    msmall[m1 >= m2] = m2[m1 >= m2]
    ilow = np.array(np.interp(m, mgrain, np.arange(nm)), dtype=int)
    jj = np.zeros((nm, nm), dtype=int) + np.arange(nm)[None, :]
    # These are the mass pairs that may lead to cancellations
    icancel = (ilow == jj)
    toobig = ilow >= nm - 1   # Sum of masses is beyond grid
    nottoobig = ilow < nm - 1    # Sum of masses is NOT beyond grid
    ilow[toobig] = nm - 2         # Handle toobig cases
    ihi = ilow + 1
    mlow = mgrain[ilow]
    mhi = mgrain[ihi]
    eps = (m - mlow) / (mhi - mlow)  # One minus Eq. A.4 of Brauer et al.
    eps[icancel] = msmall[icancel] / (mhi[icancel] - mlow[icancel])
    eps[icancel.T] = msmall[icancel.T] / (mhi[icancel.T] - mlow[icancel.T])
    eps[toobig] = 0.           # Handle toobig cases
    return ilow, eps, icancel


def create_coag_delta(mgrain, fixcancel=False):
    """
    In the Kovetz & Olund (1969) algorithm the delta(m-m1-m2) dm1 dm2 is
    numerically represented by a third-rank "tensor" that smears out over
    the two nearest mass grid points around m=m1+m2. This tensor is computed
    and returned here. It is returned in a full (non-sparse) 3D array, which
    is not particularly efficient, but easy to handle.

    See Brauer, Dullemond & Henning (2008), Appendix A1. There it is called 
    C_ijk. Note also that eps is defined oppositely to Brauer (eps --> 1-eps). 
    Here we call this tensor D_ijk (instead of Brauer's C_ijk) because it
    represents a discrete version of the Dirac delta function delta(m'+m''-m).

    Here "low" means the mass grid point just below m, and "hi" means 
    the mass grid point just above m.

    ARGUMENTS:

      mgrain       Array of mass grid points
      fixcancel    If true, then we subtract 1 from Dkij where the gain and
                   loss terms can nearly cancel each other. This avoids 
                   numerical problems if m_i<<<m_j or vv, but you need to
                   then also subtract 1 from the corresponding loss term.

    RETURNS:

      Dkij         The Kovetz & Olund (1969) delta function array.
    """
    if (fixcancel):
        assert (mgrain[1:] / mgrain[:-1]
                ).max() < 2, "Error: mass grid too coarse."
    nm = len(mgrain)
    ilow, eps, icancel = create_ilow_eps_icancel(mgrain)
    ihi = ilow + 1
    Dkij = np.zeros((nm, nm, nm))
    for i in range(nm):
        for j in range(nm):
            Dkij[ilow[i, j], i, j] = 1 - eps[i, j]
            Dkij[ihi[i, j], i, j] = eps[i, j]
            if (fixcancel):
                if icancel[i, j] or icancel[j, i]:
                    Dkij[ilow[i, j], i, j] = -eps[i, j]
    if (not fixcancel):
        icancel[:, :] = False
    return Dkij, icancel
