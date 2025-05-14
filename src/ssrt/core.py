import numpy as np

def S2RTR_DiffuseUB(eps, f, s, a, kappa_e, d, theta_deg):
    """
    Computes sigma_0_vv and sigma_0_hh for a weakly scattering Rayleigh layer 
    with albedo a < 0.2 and a PRISM-modeled ground surface.

    Parameters:
        eps : complex
            Dielectric constant of ground surface.
        f : float
            Frequency (GHz).
        s : float
            RMS height of ground surface (m).
        a : float
            Single-scattering albedo (unitless, < 0.2).
        kappa_e : float
            Extinction coefficient of Rayleigh layer (Np/m).
        d : float
            Thickness of Rayleigh layer (m).
        theta_deg : float
            Incidence angle (degrees).

    Returns:
        sigma_0_vv : float
            Backscatter coefficient for VV polarization (dB).
        sigma_0_hh : float
            Backscatter coefficient for HH polarization (dB).
    """
    theta_rad = np.radians(theta_deg)  # Convert degrees to radians
    kappa_s = a * kappa_e              # Scattering coefficient

    # --- Call PRISM-1 surface scattering model ---
    sig_s_vv_db, sig_s_hh_db, _ = PRISM1_ForwardModel(eps, theta_deg, s, f)
    sig_s_vv = 10 ** (sig_s_vv_db / 10)
    sig_s_hh = 10 ** (sig_s_hh_db / 10)

    # --- Transmissivity in the Rayleigh layer ---
    tau = kappa_e * d / np.cos(theta_rad)
    T = np.exp(-tau)

    # --- Reflectivity at the lower interface ---
    _, _, gammah, gammav, *_ = ReflTransm_PlanarBoundary(1, eps, theta_deg)

    # --- Total backscatter according to Eq. 11.23 ---
    cos_theta = np.cos(theta_rad)
    T2 = T ** 2
    sigma_0_vv = (T2 * sig_s_vv +
                  0.75 * a * cos_theta * (1 - T2) * (1 + gammav ** 2 * T2) +
                  6 * kappa_s * d * gammav * T2)

    sigma_0_hh = (T2 * sig_s_hh +
                  0.75 * a * cos_theta * (1 - T2) * (1 + gammah ** 2 * T2) +
                  6 * kappa_s * d * gammah * T2)

    # Convert to dB
    sigma_0_vv_db = 10 * np.log10(sigma_0_vv)
    sigma_0_hh_db = 10 * np.log10(sigma_0_hh)

    return sigma_0_vv_db, sigma_0_hh_db


def PRISM1_ForwardModel(eps, theta_deg, s, f):
    """
    Computes sigma_0 for all three polarization combinations 
    based on the PRISM-1 forward model.

    Parameters:
        eps : complex
            Complex dielectric constant (eps' - j*eps'').
        theta_deg : float
            Incidence angle in degrees.
        s : float
            RMS surface height (m).
        f : float
            Frequency (GHz).

    Returns:
        sig_0_vv : float
            Sigma_0 for VV polarization (dB).
        sig_0_hh : float
            Sigma_0 for HH polarization (dB).
        sig_0_hv : float
            Sigma_0 for HV polarization (dB).
    """
    theta_rad = np.radians(theta_deg)
    ks = s * (2 * np.pi * f / 0.3)  # Roughness parameter

    gamma0 = Fresn_Refl0(eps)       # Normal incidence reflectivity
    gammav, gammah = Fresn_Refl(eps, theta_rad)  # Angular-dependent reflectivity

    p = (1 - (2 * theta_rad / np.pi) ** (1 / (3 * gamma0)) * np.exp(-ks)) ** 2
    q = 0.23 * np.sqrt(gamma0) * (1 - np.exp(-ks))
    g = 0.70 * (1 - np.exp(-0.65 * ks ** 1.8))

    cos_theta = np.cos(theta_rad)
    sigvv = g * (cos_theta ** 3) / np.sqrt(p) * (gammav + gammah)

    sig_0_vv = 10 * np.log10(sigvv)
    sig_0_hh = 10 * np.log10(sigvv * p)
    sig_0_hv = 10 * np.log10(sigvv * q)

    return sig_0_vv, sig_0_hh, sig_0_hv


def Fresn_Refl0(eps):
    """
    Fresnel reflectivity at normal incidence.
    
    Parameters:
        eps : complex
            Complex dielectric constant.
    
    Returns:
        gamma0 : float
            Reflectivity (unitless).
    """
    sqrt_eps = np.sqrt(eps)
    gamma0 = np.abs((1 - sqrt_eps) / (1 + sqrt_eps)) ** 2
    return gamma0


def Fresn_Refl(eps, theta_rad):
    """
    Fresnel reflectivities for vertical and horizontal polarization.

    Parameters:
        eps : complex
            Complex dielectric constant.
        theta_rad : float
            Incidence angle in radians.

    Returns:
        gammav : float
            Vertical polarization reflectivity.
        gammah : float
            Horizontal polarization reflectivity.
    """
    rho_v, rho_h = refl_coef(theta_rad, 1, eps)
    gammav = np.abs(rho_v) ** 2
    gammah = np.abs(rho_h) ** 2
    return gammav, gammah


def refl_coef(theta_rad, eps1, eps2):
    """
    Computes vertical and horizontal polarized reflection coefficients
    of a plane dielectric surface.

    Parameters:
        theta_rad : float
            Incidence angle in radians.
        eps1 : float
            Dielectric constant of incident medium.
        eps2 : complex
            Dielectric constant of transmission medium.

    Returns:
        rho_v : complex
            Vertical reflection coefficient.
        rho_h : complex
            Horizontal reflection coefficient.
    """
    n1 = np.sqrt(eps1)
    n2 = np.sqrt(eps2)
    sin_theta = np.sin(theta_rad)
    sin_theta_ratio = (n1 * sin_theta) / n2
    costh2 = np.sqrt(1 - sin_theta_ratio ** 2)

    rho_v = -(n2 * np.cos(theta_rad) - n1 * costh2) / (n2 * np.cos(theta_rad) + n1 * costh2)
    rho_h =  (n1 * np.cos(theta_rad) - n2 * costh2) / (n1 * np.cos(theta_rad) + n2 * costh2)

    return rho_v, rho_h


def ReflTransm_PlanarBoundary(eps1, eps2, theta1d):
    """
    Computes the reflection and transmission coefficients, reflectivities, 
    and transmissivities at a planar boundary for both horizontal (h) and 
    vertical (v) polarizations.

    Parameters:
        eps1 : complex
            Relative dielectric constant of medium 1.
        eps2 : complex
            Relative dielectric constant of medium 2.
        theta1d : float
            Incidence angle in medium 1 (degrees).

    Returns:
        rhoh : complex
            Reflection coefficient for horizontal polarization.
        rhov : complex
            Reflection coefficient for vertical polarization.
        gammah : float
            Reflectivity for horizontal polarization.
        gammav : float
            Reflectivity for vertical polarization.
        tauh : complex
            Transmission coefficient for horizontal polarization.
        tauv : complex
            Transmission coefficient for vertical polarization.
        Th : float
            Transmissivity for horizontal polarization.
        Tv : float
            Transmissivity for vertical polarization.
    """
    theta1 = np.radians(theta1d)

    sqrt_eps1 = np.sqrt(eps1)
    sqrt_eps2 = np.sqrt(eps2)

    sin_theta2 = sqrt_eps1 / sqrt_eps2 * np.sin(theta1)
    cos_theta2 = np.sqrt(1 - sin_theta2 ** 2)

    rhoh = (sqrt_eps1 * np.cos(theta1) - sqrt_eps2 * cos_theta2) / \
           (sqrt_eps1 * np.cos(theta1) + sqrt_eps2 * cos_theta2)
    
    rhov = (sqrt_eps1 * cos_theta2 - sqrt_eps2 * np.cos(theta1)) / \
           (sqrt_eps1 * cos_theta2 + sqrt_eps2 * np.cos(theta1))

    tauh = 1 + rhoh
    tauv = (1 + rhov) * (np.cos(theta1) / cos_theta2)

    gammah = np.abs(rhoh) ** 2
    gammav = np.abs(rhov) ** 2

    Th = 1 - gammah
    Tv = 1 - gammav

    return rhoh, rhov, gammah, gammav, tauh, tauv, Th, Tv
