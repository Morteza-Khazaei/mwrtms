import numpy as np
from ..utils.fresnel import Fresn_Refl0, Fresn_Refl
from ..utils.util import toLambda, toPower, toDB

class PRISM1:
    """PRISM-1 rough-surface backscatter model.

    Parameters
    ----------
    f : float
        Frequency in GHz.
    theta_i : float
        Incidence angle in degrees.
    eps : complex
        Complex relative permittivity of the surface.
    s : float
        Surface RMS height (m).

    Examples
    --------
    >>> prism = PRISM1(f=5.405, theta_i=40.0, eps=4.5 + 0.5j, s=0.015)
    >>> prism.calc_sigma()  # doctest: +SKIP
    """

    def __init__(self, f, theta_i, eps, s):
        """Initialise the PRISM-1 model with scene parameters.

        Parameters
        ----------
        f : float
            Frequency in GHz.
        theta_i : float
            Incidence angle in degrees.
        eps : complex
            Complex relative permittivity of the surface.
        s : float
            Surface RMS height (m).

        Notes
        -----
        Use :meth:`calc_sigma` to obtain the four backscatter coefficients.

        """
        self.theta_rad = np.radians(theta_i)
        self.eps = eps
        lambda_m = toLambda(f)  # Wavelength in meters
        k = 2 * np.pi / lambda_m  # Wavenumber
        self.ks = k * s  # Roughness parameter
    

    def calc_sigma(self, todB=True):

        gamma0 = Fresn_Refl0(self.eps)       # Normal incidence reflectivity
        gammav, gammah = Fresn_Refl(self.eps, self.theta_rad)  # Angular-dependent reflectivity

        p = (1 - (2 * self.theta_rad / np.pi) ** (1 / (3 * gamma0)) * np.exp(-self.ks)) ** 2
        q = 0.23 * np.sqrt(gamma0) * (1 - np.exp(-self.ks))
        g = 0.70 * (1 - np.exp(-0.65 * self.ks ** 1.8))

        cos_theta = np.cos(self.theta_rad)
        sigvv = g * (cos_theta ** 3) / np.sqrt(p) * (gammav + gammah)

        sig_0_vv = sigvv
        sig_0_hh = sigvv * p
        sig_0_hv = sigvv * q
        sig_0_vh = sigvv * q

        # Convert to dB
        if todB:
            sig_0_vv = toDB(sig_0_vv)
            sig_0_hh = toDB(sig_0_hh)
            sig_0_hv = toDB(sig_0_hv)
            sig_0_vh = toDB(sig_0_vh)

        return sig_0_vv, sig_0_hh, sig_0_hv, sig_0_vh
