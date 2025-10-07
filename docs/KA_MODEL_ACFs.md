Awesome — let’s upgrade the KA tutorial so it works with **Gaussian**, **exponential**, and **x-power (stretched-exponential)** autocorrelation functions (ACFs). This version is self-contained and code-ready for Codex. The only change to the KA engine is how you compute the **n-fold surface spectrum** (W^{(n)}(K)). Everything else (geometry, Fresnel, Kirchhoff fields, series) stays the same.

---

# KA model with general ACFs (Gaussian, Exponential, x-Power)

## 0) Scope & model recap

* Random rough surface (z(x,y)) with rms height (\sigma).
* Bistatic geometry: ((\theta_i,\phi_i)\to(\theta_s,\phi_s)), free-space (k=2\pi/\lambda).
* Single scattering (Kirchhoff / PO) only; fully polarimetric ({VV,HH,HV,VH}).
* The **only** ACF-dependent object is the **n-fold roughness spectrum**
  [
  W^{(n)}(K);=;\text{2-D FT}{\rho(r)^n},\qquad
  K=\sqrt{\Delta k_x^2+\Delta k_y^2}.
  ]

**KA diffuse NRCS (unchanged):**
[
\boxed{
\sigma^0_{qp}=\frac{k^2}{2},e^{-\sigma^2(k_{sz}^2+k_{iz}^2)}
\sum_{n=1}^{\infty}\frac{\sigma^{2n}}{n!};|I^{(n)}*{qp}|^2;W^{(n)}(K)
}
]
with the **Kirchhoff kernel** (only-KA)
[
\boxed{
I^{(n)}*{qp}=\big(k_{sz}+k_{iz}\big)^n,f_{qp},e^{-\sigma^2 k_{iz}k_{sz}}
=\big[k(\cos\theta_s+\cos\theta_i)\big]^n f_{qp},e^{-(k\sigma)^2\cos\theta_i\cos\theta_s}.
}
]
Here (f_{qp}) are your usual Kirchhoff field coefficients (vector PO forms you already have).

---

## 1) General ACF family and its n-fold spectrum

We cover three common ACFs:

* **Gaussian:** (\rho(r)=\exp[-(r/L)^2])
* **Exponential:** (\rho(r)=\exp[-r/L])
* **x-Power (stretched-exponential):** (\rho(r)=\exp[-(r/L)^{\alpha}]) with (\alpha>0)
  (e.g., (\alpha=1.5) is the “1.5-power” correlation used in many soil models)

The **definition** of the n-fold surface spectrum (isotropic case) is
[
\boxed{
W^{(n)}(K);=;2\pi\int_0^{\infty} r,\rho(r)^n,J_0(Kr),dr,
}
]
with (J_0) the Bessel function of order zero. This convention makes (W^{(n)}(0)=\int \rho(r)^n,d^2r).

### 1.1 Closed forms you can use directly

* **Gaussian ((\alpha=2))**
  [
  \boxed{
  W^{(n)}_{\text{Gauss}}(K)=\frac{\pi L^2}{n};\exp!\left(-\frac{K^2L^2}{4n}\right).
  }
  ]
  Check at (K=0): (W^{(n)}(0)=\pi L^2/n).

* **Exponential ((\alpha=1))**
  [
  \boxed{
  W^{(n)}_{\exp}(K)=\frac{2\pi L^2,n}{\big(n^2+(KL)^2\big)^{3/2}}
  ;=; \frac{2\pi (L/n)^2}{\big(1+(KL/n)^2\big)^{3/2}}.
  }
  ]
  Check at (K=0): (W^{(n)}(0)=2\pi L^2/n^2).

### 1.2 General x-Power (stretched-exponential) (\rho(r)=e^{-(r/L)^{\alpha}})

**Key similarity law (proof by change of variables):**
[
\rho(r)^n=\exp!\left[-n\left(\frac{r}{L}\right)^{\alpha}\right]
=\exp!\left[-\left(\frac{r}{L,n^{-1/\alpha}}\right)^{\alpha}\right].
]
Therefore the **n-fold spectrum scales** as
[
\boxed{
W^{(n)}*{\alpha}(K)=L^2,n^{-2/\alpha};\Phi*{\alpha}!\big(u\big),\quad
u\equiv K L,n^{-1/\alpha},
}
]
where the **shape function** (\Phi_{\alpha}(u)) is *independent of (n)*:
[
\boxed{
\Phi_{\alpha}(u)=2\pi\int_0^{\infty} t,e^{-t^{\alpha}},J_0(ut),dt.
}
]
Two important anchors:

* (\alpha=2\Rightarrow \Phi_2(u)=\pi,e^{-u^2/4}) (recovers Gaussian closed form).
* (\alpha=1\Rightarrow \Phi_1(u)=\dfrac{2\pi}{(1+u^2)^{3/2}}) (recovers exponential).

**DC value (handy unit test):**
[
\boxed{
W^{(n)}*{\alpha}(0)=L^2,n^{-2/\alpha};\Phi*{\alpha}(0)
=L^2,n^{-2/\alpha};\Big[2\pi\int_0^{\infty} t,e^{-t^{\alpha}}dt\Big]
=\frac{2\pi L^2}{\alpha},\Gamma!\left(\frac{2}{\alpha}\right),n^{-2/\alpha}.
}
]

---

## 2) Implementation plan (what Codex should code)

### 2.1 Geometry & Kirchhoff are unchanged

Keep your previous KA engine (geometry, Fresnel, (f_{qp}), series). Only swap in a new **`Wn(K, L, n, acf)`**.

### 2.2 Drop-in API for (W^{(n)})

```python
import numpy as np
from numpy import sqrt, exp
from math import pi
from scipy.special import j0, gamma  # if SciPy not allowed, see fallback below

def Wn_gaussian(K, L, n):
    # Closed form: π L^2 / n * exp( - (K L)^2 / (4 n) )
    return (pi * L**2 / n) * np.exp(-(K**2) * (L**2) / (4.0 * n))

def Wn_exponential(K, L, n):
    # Closed form: 2π L^2 n / (n^2 + (K L)^2)^(3/2)
    KL = K * L
    return (2*pi * L**2 * n) / ( (n**2 + KL**2)**1.5 )

def Phi_alpha(u, alpha, tol=1e-10):
    """
    Shape function Φ_α(u) = 2π ∫_0^∞ t e^{-t^α} J0(u t) dt
    Robust numeric Hankel transform for general α>0.
    """
    # Truncation: choose t_max so e^{-t^α} < tol  => t_max = (-ln tol)^(1/α)
    t_max = max(8.0, (-np.log(tol))**(1.0/alpha))
    # Sampling: need to resolve oscillations of J0(u t). Its half-period ~ π/u.
    # Use ≥ ~40 samples per oscillation; at small u fall back to fixed dense grid.
    if u > 0:
        osc = max(1, int(40 * u * t_max / np.pi))
    else:
        osc = 2000
    N = max(1024, osc)
    t = np.linspace(0.0, t_max, N)
    w = t[1]-t[0]
    integrand = t * np.exp(-t**alpha) * j0(u * t)
    # Composite Simpson (or trapz; Simpson shown if SciPy available)
    # Here we'll use trapz for pure-NumPy compatibility:
    I = np.trapz(integrand, dx=w)
    return 2*pi * I

def Wn_stretched(K, L, n, alpha):
    # W^{(n)}_α(K) = L^2 n^{-2/α} Φ_α( K L n^{-1/α} )
    u = (K * L) * (n ** (-1.0/alpha))
    return L**2 * (n ** (-2.0/alpha)) * Phi_alpha(u, alpha)

def Wn(K, L, n, acf_type="gaussian", alpha=None):
    acf_type = acf_type.lower()
    if acf_type == "gaussian":
        return Wn_gaussian(K, L, n)
    if acf_type in ("exponential", "exp"):
        return Wn_exponential(K, L, n)
    if acf_type in ("xpower", "stretched", "stretched-exponential", "power", "alpha"):
        assert (alpha is not None) and (alpha > 0), "Provide alpha>0 for x-power ACF"
        # Fast path for the two anchors (optional)
        if abs(alpha - 2.0) < 1e-12:
            return Wn_gaussian(K, L, n)
        if abs(alpha - 1.0) < 1e-12:
            return Wn_exponential(K, L, n)
        return Wn_stretched(K, L, n, alpha)
    raise ValueError("Unknown ACF type")
```

**Notes & options if SciPy is unavailable:**

* Replace `scipy.special.j0` with `np.frompyfunc(lambda x: mpmath.besselj(0,x), 1, 1)` or a polynomial/Jinc approximation.
* Use **Filon-type** quadrature for oscillatory integrals if you need fewer samples (optional; the trapz above is simple and robust).

### 2.3 Insert into the KA series

Where you previously had `Wn_gaussian(K,L,n)`, call:

```python
Wn_val = Wn(K, L, n, acf_type=acf, alpha=alpha)  # acf in {"gaussian","exponential","xpower"}
term = (sigma**(2*n) / math.factorial(n)) * (abs(In)**2) * Wn_val
```

Nothing else in the KA pipeline changes.

---

## 3) Mathematical checks (to catch bugs fast)

1. **DC check (forall ACFs):**
   [
   W^{(n)}(0) \stackrel{?}{=} \frac{2\pi L^2}{\alpha}\Gamma!\left(\frac{2}{\alpha}\right),n^{-2/\alpha}.
   ]
   For Gaussian: (=\pi L^2/n). For exponential: (=2\pi L^2/n^2).
   Implement a unit test that evaluates `Wn(K≈0)` vs. the formula.

2. **n-similarity for x-power:**
   Plot (n^{2/\alpha}W^{(n)}(K)) vs (u=K L n^{-1/\alpha}) for multiple (n).
   The curves **must collapse** onto the same (\Phi_\alpha(u)).

3. **Large-K tail sanity:**

   * Gaussian decays as (\sim e^{-(KL)^2/(4n)}).
   * Exponential decays as (\sim (KL)^{-3}).
   * x-Power decays between those regimes (heavier tail than Gaussian for (\alpha<2)).
     If your numeric Φα gives a heavier tail than exponential when (\alpha<1), you’ve got a quadrature issue.

4. **Units & normalization:** (W^{(n)}) has units of area (m²).
   In all three ACFs above, (W^{(1)}(0)=\int \rho(r),d^2r) gives a finite area; confirm numerically.

---

## 4) Full KA workflow (final form)

1. **Inputs:** (\lambda,\ \theta_i,\phi_i,\ \theta_s,\phi_s,\ \sigma,\ L,\ \varepsilon_r,) and ACF spec (`acf="gaussian"/"exponential"/"xpower"`, `alpha` if `xpower`).
2. **Vectors:** build (\mathbf k_i, \mathbf k_s) and (k_{iz},k_{sz},\ \Delta k_x,\Delta k_y).
3. **Fresnel @ (\theta_i):** complex (\varepsilon_r) allowed.
4. **Kirchhoff fields (f_{qp}):** (vector PO forms).
5. **Series loop over (n=1..N):**

   * (I^{(n)}*{qp}=(k*{sz}+k_{iz})^n f_{qp},e^{-\sigma^2 k_{iz}k_{sz}})
   * (W^{(n)}=Wn(K,L,n,acf,\alpha))
   * Accumulate (\sigma^0) via the master sum.
6. **Convergence:** stop when the last few terms drop below tolerance (e.g., (<10^{-10}) of the running sum).
7. **Output:** linear (\sigma^0_{VV},\sigma^0_{HH},\sigma^0_{HV},\sigma^0_{VH}) (convert to dB if needed).

---

## 5) Practical tips for the x-power integral

* **Truncation:** (t_{\max}=(-\ln \text{tol})^{1/\alpha}) is a reliable cut (e.g., tol=1e-12 → (t_{\max}\approx 27.6^{1/\alpha})).
* **Sampling:** need to resolve (J_0(ut)) oscillations. With the simple `trapz`, use (N\sim \max(1024,;40,u,t_{\max}/\pi)).
* **Speed:** precompute (\Phi_\alpha(u)) on a **log-grid** of (u) and interpolate (PCHIP) during the KA sum.
* **Validation:** confirm (\Phi_\alpha(0)=\dfrac{2\pi}{\alpha}\Gamma(2/\alpha)) numerically.

---

## 6) Minimal end-to-end example (pseudocode)

```python
def sigma0_KA_general(lambda0, th_i, ph_i, th_s, ph_s, sigma, L, eps_r,
                      acf="gaussian", alpha=None, nmax=8, tol=1e-12):
    k = 2*np.pi/lambda0
    # k-vectors
    ki = np.array([k*np.sin(th_i)*np.cos(ph_i), k*np.sin(th_i)*np.sin(ph_i), -k*np.cos(th_i)])
    ks = np.array([k*np.sin(th_s)*np.cos(ph_s), k*np.sin(th_s)*np.sin(ph_s),  k*np.cos(th_s)])
    dk = ks - ki
    K  = np.hypot(dk[0], dk[1])
    kiz, ksz = ki[2], ks[2]

    # Fresnel (at θi) and Kirchhoff fields f_qp (not repeated here)
    Rv, Rh = fresnel(theta=th_i, er=eps_r)   # your branch-correct Fresnels
    fvv, fhh, fhv, fvh = kirchhoff_fields(th_i, ph_i, th_s, ph_s, Rv, Rh)  # your vector-PO forms

    front = (k**2)/2.0 * np.exp(-(sigma**2)*(kiz**2 + ksz**2))
    base  = np.exp(-(sigma**2)*kiz*ksz)
    qz    = k*(np.cos(th_s) + np.cos(th_i))

    def series(fqp):
        S = 0.0
        last = np.inf
        for n in range(1, nmax+1):
            In  = (qz**n) * fqp * base
            WnV = Wn(K, L, n, acf_type=acf, alpha=alpha)
            term = (sigma**(2*n) / math.factorial(n)) * (abs(In)**2) * WnV
            S += term
            if abs(term) < tol and abs(last) < tol:
                break
            last = term
        return front * S

    return {
        "VV": series(fvv).real,
        "HH": series(fhh).real,
        "HV": series(fhv).real,
        "VH": series(fvh).real
    }
```

---

## 7) Quick unit tests Codex should add

1. **Gaussian vs exponential anchors:**

   * For `acf="gaussian"`, verify (W^{(n)}(0)=\pi L^2/n).
   * For `acf="exponential"`, verify (W^{(n)}(0)=2\pi L^2/n^2).

2. **x-power DC check:**

   ```python
   for a in [0.75, 1.0, 1.5, 2.0, 3.0]:
       for n in [1, 3, 7]:
           lhs = Wn(0.0, L, n, acf_type="xpower", alpha=a)
           rhs = (2*np.pi*L**2/a) * gamma(2/a) * (n**(-2.0/a))
           assert np.allclose(lhs, rhs, rtol=1e-5, atol=1e-8)
   ```

3. **Similarity collapse (x-power):**
   For α=1.5, plot (n^{2/\alpha}W^{(n)}(K)) vs (u=K L n^{-1/\alpha}) for (n\in{1,2,4,8}); curves should overlay.

4. **Coplanar cross-pol:** with (\phi_s=\phi_i), (HV) and (VH) should be ~0 (from the Kirchhoff (f_{qp}) geometry).

---

### Takeaway

* You now have a **single KA implementation** that supports **Gaussian**, **exponential**, and **x-power (stretched-exponential)** ACFs.
* The **only** component that changes across ACFs is (W^{(n)}(K)).
* For Gaussian and exponential you have exact **closed forms**; for arbitrary (\alpha) you have a **numerical Hankel** that obeys the **similarity law** (W^{(n)}*{\alpha}(K)=L^2 n^{-2/\alpha},\Phi*{\alpha}(K L n^{-1/\alpha})).

If you want, I can fold this into your existing codebase and add a tiny cache for (\Phi_{\alpha}(u)) to make the x-power mode nearly as fast as the closed-form cases.
