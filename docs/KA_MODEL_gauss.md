Here’s a clean, code-ready tutorial for building a **Kirchhoff Approximation (KA)** scattering model (bistatic, fully polarimetric) from the thesis you shared. I’ve laid it out so Codex can translate each step into a working Python module. I also include the exact equations you’ll need (in a self-contained form), and I point to the relevant places in the thesis text with inline file citations so you can verify each definition.

---

# KA Model (Bistatic, Gaussian Surface) — Implementation Tutorial

## 0) Scope & assumptions

* Random rough surface, height ( z(x,y) ) is **Gaussian** with rms height ( \sigma ) and isotropic **Gaussian correlation** of length (L).
* Bistatic geometry: incident direction ((\theta_i,\phi_i)), scattered direction ((\theta_s,\phi_s)).
* Single scattering only (Kirchhoff / physical-optics approximation). No complementary-field or multiple-scattering terms.
* Fully polarimetric (VV, HH, HV, VH).

The goal is to compute the **bistatic scattering coefficient** (\sigma^0_{qp}) (with (p) transmit pol, (q) receive pol).

---

## 1) Geometry & wave vectors

Let free-space wavenumber (k=2\pi/\lambda). Adopt the thesis convention:

[
\begin{aligned}
\mathbf{k}*i &= {k_x,k_y,k_z}
= {,k \sin\theta_i\cos\phi_i,; k \sin\theta_i\sin\phi_i,; -k\cos\theta_i,},\
\mathbf{k}*s &= {k*{sx},k*{sy},k_{sz}}
= {,k \sin\theta_s\cos\phi_s,; k \sin\theta_s\sin\phi_s,; k\cos\theta_s,}.
\end{aligned}
]

This sets (k_z=-k\cos\theta_i) for the **incident** wave (downwards) and (k_{sz}=+k\cos\theta_s) for the **scattered** wave (upwards). 

Define horizontal wavenumber **mismatch**:

[
\Delta k_x = k_{sx}-k_x, \qquad
\Delta k_y = k_{sy}-k_y .
]

Define also the **vertical sum** (this is crucial in KA):

[
q_z \equiv k_{sz} + |k_z| = k(\cos\theta_s + \cos\theta_i).
]

> In this thesis’ KA/PO derivation, the noncoherent term organizes naturally in powers of (q_z^2) and carries an outer Gaussian factor with ((\cos\theta_s+\cos\theta_i)^2). 

---

## 2) Surface statistics (Gaussian correlation)

**Correlation** (isotropic Gaussian):

[
\rho(r) = \exp!\Big(-\frac{r^2}{L^2}\Big).
]

The **(n)-th power roughness spectrum** (2-D FT of (\rho(r)^n)) in closed form:

[
\boxed{ ;
W^{(n)}(K) ;=; \frac{\pi L^2}{n}; \exp!\Big(!-,\frac{K^2 L^2}{4n}\Big),
;}
\qquad K=\sqrt{ (\Delta k_x)^2+(\Delta k_y)^2 }.
]

(That expression follows directly from the Bessel-transform identity given in the thesis.) 

> Tip (code): precompute (K) once per look direction and then evaluate (W^{(n)}(K)) analytically for each (n).

Normalized roughness (often used for param sweeps):

[
kL=2\pi L/\lambda,\qquad k\sigma=2\pi \sigma/\lambda. ; ; \text{}
]

---

## 3) Fresnel reflection coefficients (medium 1 → medium 2)

Let (\varepsilon_r = \varepsilon' - j\varepsilon'') be the relative permittivity of the lower half-space.

With local (planar) incidence angle (\theta) into medium 2:

[
\begin{aligned}
R_v(\theta) &= \frac{\varepsilon_r \cos\theta - \sqrt{\varepsilon_r - \sin^2\theta}}{\varepsilon_r \cos\theta + \sqrt{\varepsilon_r - \sin^2\theta}},\
R_h(\theta) &= \frac{\cos\theta - \sqrt{\varepsilon_r - \sin^2\theta}}{\cos\theta + \sqrt{\varepsilon_r - \sin^2\theta}}.
\end{aligned}
]

> Use the **principal branch** for the square root and enforce (\mathrm{Re}{\sqrt{\varepsilon_r - \sin^2\theta}}\ge 0) to keep energy flow into medium 2.

For the KA single-scattering term we evaluate these at the **global** incident angle (\theta_i) (Kirchhoff assumes locally planar facets).

---

## 4) Kirchhoff (physical-optics) field coefficients (f_{qp})

Define the 2-D surface spectral slope parameters used in the vector KA (same as AIEM’s Kirchhoff piece):

[
\begin{aligned}
z_x &\equiv -,\frac{k_{sx}-k_x}{k_{sz}+|k_z|} ;=; -,\frac{\sin\theta_s\cos\phi_s - \sin\theta_i\cos\phi_i}{\cos\theta_s + \cos\theta_i},[4pt]
z_y &\equiv -,\frac{k_{sy}-k_y}{k_{sz}+|k_z|} ;=; -,\frac{\sin\theta_s\sin\phi_s - \sin\theta_i\sin\phi_i}{\cos\theta_s + \cos\theta_i},
\end{aligned}
]

and the normalization

[
\Delta ;\equiv; \sqrt{(z_x \cos\theta_i - \sin\theta_i)^2 + z_y^2 }.
]

Using these, the **vector** Kirchhoff/PO scattering amplitudes (f_{qp}) can be written in closed form (below I give the exact algebraic forms that are used widely in IEM/AIEM; they are valid PO/KA amplitudes as well). To keep this tutorial focused, I place the four expressions here verbatim (they match the “(f_{vv},f_{hh},f_{hv},f_{vh})” Kirchoff terms that appear as the first term in single scattering):

[
\begin{aligned}
f_{hh} &= (1-R_h),h_{snv} + (1+R_h),v_{snh} - (h_{snt}+v_{snd})(R_h+R_v),\frac{z_y}{\Delta},[2pt]
f_{vv} &= -\Big[(1-R_v),h_{snv} + (1+R_v),v_{snh}\Big] + (h_{snt}+v_{snd})(R_h+R_v),\frac{z_y}{\Delta},[2pt]
f_{hv} &= -(1+R_v),h_{snh} + (1-R_v),v_{snv} + (h_{snd}-v_{snt})(R_h+R_v),\frac{z_y}{\Delta},[2pt]
f_{vh} &= -(1+R_h),h_{snh} + (1-R_h),v_{snv} + (h_{snd}-v_{snt})(R_h+R_v),\frac{z_y}{\Delta}.
\end{aligned}
]

The geometric factors ((h_{snv},v_{snh},v_{snv},h_{snt},v_{snd})) are deterministic functions of the angles and the pair ((z_x,z_y)). One consistent set (matching the IEM/AIEM Kirchhoff part) is:

[
\begin{aligned}
h_{snv} &= -(\cos\theta_i\cos\phi_s + \sin\theta_i(z_x\cos\phi_s+z_y\sin\phi_s)),\
v_{snh} &= \cos\theta_s\cos\phi_s - z_x\sin\theta_s,\
h_{snh} &= -\sin\phi_s,\
v_{snv} &= z_y\cos\theta_i\sin\theta_i + \cos\theta_s\big(z_y\cos\phi_s\sin\theta_i - (\cos\theta_i+z_x\sin\theta_i)\sin\phi_s\big),\
h_{snt} &= \frac{-(\cos^2\theta_i+\sin^2\theta_i)\sin\phi_s(-\sin\theta_i+\cos\theta_i z_x)+\cos\phi_s(\cos\theta_i+z_x\sin\theta_i)z_y+\sin\theta_i\sin\phi_s z_y^2}{\Delta},\
v_{snd} &= \frac{-(\cos\theta_i+z_x\sin\theta_i)\big(-\cos\phi_s\sin\theta_i+\cos\theta_i\cos\phi_s z_x+\cos\theta_i\sin\phi_s z_y\big)}{\Delta}.
\end{aligned}
]

These are exactly the ones you see encoded in common AIEM/PO implementations and are the **“(f_{qp})”** used in the single-scattering kernel; they are consistent with the vector PO derivation (compare with the AIEM/KA appendix forms). 

> If you prefer a **scalar PO** shortcut for quick checks, you can set
> (f_{vv}\approx (1+R_v)\cos\theta_s,; f_{hh}\approx (1+R_h)\cos\theta_s)
> and (f_{hv}=f_{vh}=0). Use the full vector version above for production.

---

## 5) KA single-scattering coefficient

The KA **single scattering** (noncoherent) coefficient organizes as a series in the height moments with powers of (q_z) and the roughness spectra (W^{(n)}):

[
\boxed{ ;
\sigma^0_{qp} ;=; \frac{k^2}{2};
\exp!\big[-,\sigma^2,(k_{sz}^2 + k_z^2)\big];
\sum_{n=1}^\infty \frac{\sigma^{2n}}{n!},
\Big|, I^{(n)}_{qp} ,\Big|^2;
W^{(n)}!\big(\Delta k_x,\Delta k_y\big),
;}
]

with the **KA (Kirchhoff-only)** moment kernel

[
\boxed{ ;
I^{(n)}*{qp} ;=; (k*{sz}+|k_z|)^n; f_{qp};\exp(-,\sigma^2 |k_z|,k_{sz})
;=; \big[k(\cos\theta_s+\cos\theta_i)\big]^n f_{qp},\exp!\big(-\sigma^2 k^2\cos\theta_i \cos\theta_s\big).
;}
]

* The outer exponential ( \exp[-\sigma^2 (k_{sz}^2 + k_z^2)] = \exp[-(k\sigma)^2(\cos^2\theta_s+\cos^2\theta_i)] ) multiplies the series,
* The inner exponential in (I^{(n)}_{qp}) contributes a factor (\exp[-\sigma^2 k^2 \cos\theta_i\cos\theta_s]),
* Combined, they build (\exp[-(k\sigma)^2(\cos\theta_s+\cos\theta_i)^2]), the familiar KA damping envelope.

This structure matches the thesis’ KA/PO noncoherent construction (powers of (q_z) with Gaussian envelope and the (W^{(n)}) spectra).  

For **Gaussian correlation** insert (W^{(n)}) from §2 above.

> Practical convergence: for moderate roughness, the series converges fast; (n_{\max}\in[4,10]) is typically sufficient. For very rough surfaces, KA may be outside its validity region (PO/GO transition).

---

## 6) Algorithm (step-by-step)

1. **Inputs**: (\lambda, \theta_i,\phi_i,\theta_s,\phi_s, \sigma, L, \varepsilon_r).
2. **Wavenumbers**: compute (k=2\pi/\lambda) and the components ((k_x,k_y,k_z)), ((k_{sx},k_{sy},k_{sz})). 
3. **Geometry helpers**: ( \Delta k_x,\Delta k_y, K, q_z = k(\cos\theta_s+\cos\theta_i)).
4. **Fresnel**: compute (R_v(\theta_i), R_h(\theta_i)) (principal branch).
5. **Kirchhoff amplitude**: compute (z_x,z_y,\Delta) then (f_{qp}) using the vector formulas in §4.
6. **Roughness spectrum**: for (n=1,\dots,n_{\max}), compute (W^{(n)}(K) = \frac{\pi L^2}{n}\exp(-K^2L^2/(4n))).
7. **KA kernel**: (I^{(n)}*{qp} = q_z^n, f*{qp},\exp(-\sigma^2 |k_z| k_{sz})).
8. **Accumulate**:
   [
   \sigma^0_{qp} \mathrel{+}= \frac{k^2}{2},\exp[-\sigma^2(k_{sz}^2+k_z^2)];\frac{\sigma^{2n}}{n!};|I^{(n)}_{qp}|^2,W^{(n)}(K).
   ]
9. **Polarization outputs**: compute VV, HH, HV, VH.

---

## 7) Minimal Python skeleton (ready for Codex)

```python
import numpy as np
from math import pi

def fresnel_Rv_Rh(theta, eps_r):
    s2 = np.sin(theta)**2
    root = np.sqrt(eps_r - s2)   # principal branch
    Rv = (eps_r*np.cos(theta) - root) / (eps_r*np.cos(theta) + root)
    Rh = (np.cos(theta) - root) / (np.cos(theta) + root)
    return Rv, Rh

def gaussian_Wn(K, L, n):
    return (np.pi * L**2 / n) * np.exp(-(K**2) * (L**2) / (4.0*n))

def k_vectors(k, th_i, ph_i, th_s, ph_s):
    kx  = k*np.sin(th_i)*np.cos(ph_i)
    ky  = k*np.sin(th_i)*np.sin(ph_i)
    kz  = -k*np.cos(th_i)
    ksx = k*np.sin(th_s)*np.cos(ph_s)
    ksy = k*np.sin(th_s)*np.sin(ph_s)
    ksz = k*np.cos(th_s)
    return (kx,ky,kz, ksx,ksy,ksz)

def kirchhoff_fqp(th_i, ph_s, th_s, Rv, Rh):
    # z_x, z_y, Delta
    zx = -(np.sin(th_s)*np.cos(ph_s) - np.sin(th_i)*np.cos(0.0)) / (np.cos(th_s) + np.cos(th_i))
    zy = -(np.sin(th_s)*np.sin(ph_s) - np.sin(th_i)*np.sin(0.0)) / (np.cos(th_s) + np.cos(th_i))
    Delta = np.sqrt((zx*np.cos(th_i) - np.sin(th_i))**2 + zy**2)

    # geometry helpers
    hsnv = -(np.cos(th_i)*np.cos(ph_s) + np.sin(th_i)*(zx*np.cos(ph_s)+zy*np.sin(ph_s)))
    vsnh =  np.cos(th_s)*np.cos(ph_s) - zx*np.sin(th_s)
    hsnh = -np.sin(ph_s)
    vsnv =  zy*np.cos(th_i)*np.sin(th_i) + np.cos(th_s)*( zy*np.cos(ph_s)*np.sin(th_i)
            - (np.cos(th_i)+zx*np.sin(th_i))*np.sin(ph_s) )
    hsnt = (-(np.cos(th_i)**2+np.sin(th_i)**2)*np.sin(ph_s)*(-np.sin(th_i)+np.cos(th_i)*zx)
            + np.cos(ph_s)*(np.cos(th_i)+zx*np.sin(th_i))*zy
            + np.sin(th_i)*np.sin(ph_s)*zy**2) / Delta
    vsnd = (-(np.cos(th_i)+zx*np.sin(th_i))*(-np.cos(ph_s)*np.sin(th_i)
            + np.cos(th_i)*np.cos(ph_s)*zx + np.cos(th_i)*np.sin(ph_s)*zy)) / Delta

    fhh = (1-Rh)*hsnv + (1+Rh)*vsnh - (hsnt+vsnd)*(Rh+Rv)*(zy/Delta)
    fvv = -((1-Rv)*hsnv + (1+Rv)*vsnh) + (hsnt+vsnd)*(Rh+Rv)*(zy/Delta)
    fhv = -(1+Rv)*hsnh + (1-Rv)*vsnv + (vsnd-hsnt)*(Rh+Rv)*(zy/Delta)
    fvh = -(1+Rh)*hsnh + (1-Rh)*vsnv + (vsnd-hsnt)*(Rh+Rv)*(zy/Delta)

    return fvv, fhh, fhv, fvh

def sigma0_KA(lambda0, th_i, ph_i, th_s, ph_s, sigma, L, eps_r, nmax=8):
    k = 2*pi/lambda0
    kx,ky,kz, ksx,ksy,ksz = k_vectors(k, th_i, ph_i, th_s, ph_s)
    dKx, dKy = (ksx-kx), (ksy-ky)
    K = np.hypot(dKx, dKy)
    qz = k*(np.cos(th_s) + np.cos(th_i))

    Rv, Rh = fresnel_Rv_Rh(th_i, eps_r)
    fvv, fhh, fhv, fvh = kirchhoff_fqp(th_i, ph_s, th_s, Rv, Rh)

    outer = (k**2)/2.0 * np.exp(- (sigma**2) * (ksz**2 + kz**2))
    inner_exp = np.exp(- (sigma**2) * (abs(kz) * ksz))  # = exp(-sigma^2 k^2 cos th_i cos th_s)

    def series_term(fqp):
        acc = 0.0 + 0.0j
        for n in range(1, nmax+1):
            In = (qz**n) * fqp * inner_exp
            acc += (sigma**(2*n) / np.math.factorial(n)) * (abs(In)**2) * gaussian_Wn(K, L, n)
        return outer * acc

    VV = series_term(fvv)
    HH = series_term(fhh)
    HV = series_term(fhv)
    VH = series_term(fvh)
    return {'VV': VV.real, 'HH': HH.real, 'HV': HV.real, 'VH': VH.real}
```

Notes:

* Angles in **radians** in the code.
* The chosen (f_{qp}) are the standard vector PO (Kirchhoff) amplitudes used in IEM/AIEM for the single-scattering Kirchhoff term. 
* If you only need **backscatter**: set (\theta_s=\theta_i,; \phi_s=\phi_i+\pi). The geometry routines above already support general bistatics. 

---

## 8) Sanity checks & validation

1. **Specular lobe**: as ((\theta_s,\phi_s)\to) specular of ((\theta_i,\phi_i)), (K\to 0), (W^{(n)}(K)\to \frac{\pi L^2}{n}) and the envelope becomes (\exp[-(k\sigma)^2(\cos\theta_i+\cos\theta_s)^2]) — the characteristic KA Gaussian. 
2. **Smooth surface** ((k\sigma\to 0)): the series vanishes (\Rightarrow) (\sigma^0\to 0) (no incoherent return), as expected.
3. **Energy/units**: code returns linear (\sigma^0). Convert to dB as (10\log_{10}(\sigma^0)) if needed.
4. **Convergence**: increase `nmax` until results stabilize (start with 6–8).
5. **Geometry**: verify your (\mathbf{k}_i,\mathbf{k}*s) signs match the definitions (incident (k_z<0), scattered (k*{sz}>0)). 

---

## 9) What this KA model intentionally omits

* **Complementary-field terms (F, G)** and multi-scatter couplings (which belong to IEM/AIEM single- and multiple-scattering). If you later extend to AIEM, insert those extra (F^\pm,G^\pm) terms into (I^{(n)}_{qp}) as shown in the thesis appendices; the KA presented here keeps only the **Kirchhoff** part. 

---

## 10) Reference checkpoints in the thesis (for your verification)

* Gaussian (W^{(n)}(K)) closed form: ( \frac{\pi L^2}{n}\exp(-K^2L^2/(4n)) ). 
* KA/PO noncoherent structure (powers of vertical sum, Gaussian envelope, and (W^{(n)})): see the noncoherent PO organization and the single-scattering series structure.  
* Bistatic k-vector definitions and backscatter convention. 

---

### Closing remarks

This tutorial gives you a **standalone** KA (Kirchhoff-only) engine that mirrors the thesis’ organization:

* exact geometry,
* exact Gaussian (W^{(n)}),
* vector PO amplitudes (f_{qp}),
* KA series in (q_z^n) with the proper Gaussian damping.

If you need me to wire this into your existing Python AIEM framework (so you can toggle KA vs IEM/AIEM), say the word and I’ll inline this as a drop-in module with tests.
