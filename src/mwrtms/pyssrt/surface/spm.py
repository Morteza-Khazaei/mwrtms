import numpy as np

from ..utils.util import toLambda, toDB


class SPM3D:
    """First-order Small Perturbation Method (SPM) backscatter model.

    Parameters
    ----------
    fr : float
        Frequency in GHz.
    sig : float
        Surface RMS height (m).
    L : float
        Correlation length (m).
    thi : float
        Incidence angle in degrees.
    eps : complex
        Complex relative permittivity of the surface.

    Examples
    --------
    >>> spm = SPM3D(fr=5.405, sig=0.01, L=0.1, thi=40.0, eps=5.0 + 1.0j)
    >>> spm.calc_sigma()  # doctest: +SKIP
    """

    def __init__(self, fr, sig, L, thi, eps):
        """Initialise the SPM3D model with geometric and dielectric parameters."""
        self.theta_rad = np.radians(thi)
        self.eps = eps
        self.h = sig
        self.corr_length = L
        lambda_m = toLambda(fr)
        self.k = 2 * np.pi / lambda_m

    def calc_sigma(self, todB=True):
        theta = self.theta_rad
        h = self.h
        correlation_length = self.corr_length
        k = self.k
        eps = self.eps
        k1 = k * np.sqrt(eps)

        kzi = k * np.cos(theta)
        kxi = k * np.sin(theta)
        k1zi = np.sqrt(k1 ** 2 - kxi ** 2)

        kz = k * np.cos(theta)
        kx = -k * np.sin(theta)
        k1z = np.sqrt(k1 ** 2 - kx ** 2)

        termh = np.abs((kzi - k1zi) / (kzi + k1zi)) ** 2

        numerator_v = (k1 ** 2 - k ** 2) * (k1 ** 2 * k ** 2 * np.sin(theta) ** 2 + k ** 2 * k1z * k1z)
        denominator_v = (k1 ** 2 * kz + k ** 2 * k1z) ** 2
        termv = np.abs(numerator_v / denominator_v) ** 2

        krms = k * h
        sin_theta_sq = np.sin(theta) ** 2
        cos_theta = np.cos(theta)
        denom = (4 * k ** 2 * correlation_length ** 2 * sin_theta_sq + 1) ** (3 / 2)
        prefactor = 8 * k ** 4 * h ** 2 * correlation_length ** 2 / denom * cos_theta ** 4

        sig_0_hh = prefactor * termh
        sig_0_vv = prefactor * termv

        sqrt_eps_real = np.sqrt(np.maximum(np.real(eps), 0.0))
        ratio_term = (sqrt_eps_real - 1) / (sqrt_eps_real + 1)
        ratio_term = np.clip(ratio_term, a_min=0.0, a_max=None)
        hv_coeff = 0.23 * np.sqrt(ratio_term)
        sig_0_hv = hv_coeff * (1 - np.exp(-krms)) * sig_0_vv
        sig_0_vh = sig_0_hv

        min_val = np.finfo(float).tiny
        sig_0_vv = np.maximum(sig_0_vv, min_val)
        sig_0_hh = np.maximum(sig_0_hh, min_val)
        sig_0_hv = np.maximum(sig_0_hv, min_val)
        sig_0_vh = np.maximum(sig_0_vh, min_val)

        if todB:
            sig_0_vv = toDB(sig_0_vv)
            sig_0_hh = toDB(sig_0_hh)
            sig_0_hv = toDB(sig_0_hv)
            sig_0_vh = toDB(sig_0_vh)

        return sig_0_vv, sig_0_hh, sig_0_hv, sig_0_vh
