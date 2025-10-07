# Editor-in-Chief Report

**Subject:** Final, standalone specification for auditing and fixing the MATLAB AIEM, ready for a Python port (no external references required)

This report merges all reviewer findings into one clear, implementation-grade document. It is **self-contained**: every required equation is written out and every disputed point is resolved here. Hand this to Codex (or any engineer) to debug/port the model without opening any papers.

---

## 0) Symbols, coordinates, and shorthands (use verbatim)

* Free-space wavenumber: (k=2\pi/\lambda).
* RMS height: (\sigma), correlation length: (L).
* Dimensionless: (ks=k\sigma), (kl=kL).
* Relative permittivity: (\varepsilon_r=\varepsilon'_r + j,\varepsilon''_r) (complex allowed).
* Nonmagnetic soil: (\mu_r=1).

### Directions (unit vectors; surface normal is (+\hat z))

* Incident:
  [
  \hat{\mathbf k}_i=(\sin\theta_i\cos\phi_i,;\sin\theta_i\sin\phi_i,;-\cos\theta_i).
  ]
* Scattered:
  [
  \hat{\mathbf k}_s=(\sin\theta_s\cos\phi_s,;\sin\theta_s\sin\phi_s,;+\cos\theta_s).
  ]

### MATLAB-style shorthand

[
si=\sin\theta_i,\quad cs=\cos\theta_i,\quad
sis=\sin\theta_s,\quad css=\cos\theta_s,\quad
csfs=\cos\phi_s,\quad sfs=\sin\phi_s.
]

### Transverse spectral shift (drives the roughness spectrum)

[
\Delta k_x=k\big(sis\cdot csfs - si\cdot \cos\phi_i\big),\qquad
\Delta k_y=k\big(sis\cdot sfs  - si\cdot \sin\phi_i\big),
]
[
K=\sqrt{\Delta k_x^2+\Delta k_y^2},\qquad
K' \equiv KL = kl;\sqrt{(\Delta k_x/k)^2+(\Delta k_y/k)^2}.
]

> With (\phi_i=0) (as in the MATLAB), your existing (K) is correct.

---

## 1) Specular half-angle used in the transition/blending (✅ **correct in your code**)

The “local specular incidence” (\theta_{\mathrm{sp}}) is the **half of the angle between the two outward-pointing rays**: the reversed incident ray and the scattered ray. Define
[
\mathbf u=-\hat{\mathbf k}_i=(\sin\theta_i\cos\phi_i,;\sin\theta_i\sin\phi_i,;+\cos\theta_i),
\quad
\mathbf v=\hat{\mathbf k}*s=(\sin\theta_s\cos\phi_s,;\sin\theta_s\sin\phi_s,;+\cos\theta_s).
]
Then
[
\cos\psi'=\mathbf u!\cdot!\mathbf v
=\sin\theta_i\sin\theta_s\cos(\phi_s-\phi_i);+;\cos\theta_i\cos\theta_s,
]
[
\boxed{;\cos\theta*{\mathrm{sp}}=\cos\frac{\psi'}{2}
=\sqrt{\dfrac{1+\cos\psi'}{2}}
=\sqrt{\dfrac{1 + cs\cdot css - si\cdot sis\cdot \cos(\phi_s-\phi_i)}{2}}; .}
]
With (\phi_i=0):
[
\boxed{csl=\sqrt{\dfrac{1 + cs\cdot css - si\cdot sis\cdot csfs}{2}},\qquad
sil=\sqrt{1-csl^2}.}
]

**Sanity checks.**
Monostatic ((\theta_s=\theta_i,\ \phi_s=\phi_i+\pi\Rightarrow\cos(\phi_s-\phi_i)=-1)): (\cos\theta_{\mathrm{sp}}=1\Rightarrow \theta_{\mathrm{sp}}=0^\circ) (normal incidence on the specular facet).
Forward direction ((\theta_s=\theta_i,\ \phi_s=\phi_i)): (\cos\theta_{\mathrm{sp}}=|\cos\theta_i|).

> **Resolved point:** earlier a reviewer flagged a sign error. That was based on bisecting (\hat{\mathbf k}_i) with (\hat{\mathbf k}_s) (one inward, one outward). The correct bisector is between (-\hat{\mathbf k}_i) and (\hat{\mathbf k}_s). Your code already does the **correct** formula.

---

## 2) Fresnel reflection coefficients (fix the complex branch & normal incidence)

### 2.1 Transmitted vertical component (normalized by (k))

[
\boxed{;\frac{k_{tz}}{k}=\sqrt{\varepsilon_r-\sin^2\theta};,}
]
with the **branch chosen** so the transmitted field **decays into the lower half-space**:
[
\boxed{;\Im{k_{tz}}\ge 0\quad\text{(and typically }\Re{k_{tz}}\ge 0\text{)} .}
]
**Implementation tip:** compute (w=\sqrt{\varepsilon_r-\sin^2\theta}) (principal square root). If (\Im{w}<0), set (w:=-w). Use (k_{tz}/k=w).

### 2.2 Fresnel amplitudes (nonmagnetic soil)

[
\boxed{R_h(\theta)=\frac{cs - (k_{tz}/k)}{cs + (k_{tz}/k)},,\qquad
R_v(\theta)=\frac{\varepsilon_r,cs - (k_{tz}/k)}{\varepsilon_r,cs + (k_{tz}/k)}.}
]

### 2.3 Normal incidence (use the same (r_0) for both polarizations)

[
\boxed{R_h(0)=R_v(0)=\frac{1-\sqrt{\varepsilon_r}}{1+\sqrt{\varepsilon_r}};\equiv r_0.}
]
**Bug in MATLAB:** `rv0 =(sqrt(er)-1)/(sqrt(er)+1)`, `rh0 = -(sqrt(er)-1)/(sqrt(er)+1)`.
**Fix:** set **both** to (r_0=(1-\sqrt{\varepsilon_r})/(1+\sqrt{\varepsilon_r})).

---

## 3) Transition function (replace legacy block)

Blend the Fresnel coefficients at the incident angle and at the specular half-angle using
[
\boxed{R_p^{(T)};=;R_p(\theta_i);+;\Big(R_p(\theta_{\mathrm{sp}})-R_p(\theta_i)\Big);\gamma_p,\qquad p\in{h,v}.}
]
with
[
\boxed{;\gamma_p=1-\dfrac{S_p}{S_p^{(0)}};.}
]

**How to compute (S_p) and (S_p^{(0)}) (drop-in recipe):**

1. Temporarily **freeze all Fresnels** to the normal-incidence value (r_0) from §2.3: use (R_h(\cdot)=R_v(\cdot)=r_0) everywhere in the single-scattering computation.
2. Using your existing assembly (Kirchhoff + cross + complementary):

   * Compute the **total** co-pol backscatter (\sigma^{\text{all}}_{pp}).
   * Compute the **complementary-only** part (\sigma^{\text{comp}}*{pp}) (i.e., rebuild (I^{(n)}*{pp}) **without** the Kirchhoff term; keep only the eight complementary terms).
3. Define (S_p:=\sigma^{\text{comp}}_{pp}). Obtain (S_p^{(0)}) numerically by repeating step 2 at a tiny (ks) (e.g., (ks=10^{-6})); the value stabilizes rapidly.
4. Form (\gamma_p=1-S_p/S_p^{(0)}) and then (R_p^{(T)}) by the boxed formula.

**What to delete/repair in MATLAB:**

* Remove the ad-hoc (F_{t*})/`st*` block entirely.
* If a temporary scaffold is kept while migrating:

  * **Fix the typo:** in the H-path use `rh0` (not `rv0`).
  * Make **both** branches use the **same** Gaussian factor, e.g. multiply by (\exp{-(ks\cdot cs)^2}); do **not** divide by (\exp{+(ks\cdot cs)^2}) anywhere.

---

## 4) Roughness spectra (W^{(n)}(K)) (2-D FT of (\rho(r)^n))

Let (K'=KL). The model needs the **n-fold** spectrum (W^{(n)}), i.e. the 2-D Fourier transform of (\rho(r)^n).

### 4.1 Gaussian correlation (\rho(r)=\exp[-(r/L)^2])  (✅ keep)

[
\boxed{,W^{(n)}(K);\propto;\frac{L^2}{n};\exp!\left(-\frac{(K L)^2}{4n}\right)
=\frac{L^2}{n};\exp!\left(-\frac{K'^2}{4n}\right).}
]
Your implementation matches this (shape and (1/n) scaling).

### 4.2 Exponential correlation (\rho(r)=\exp[-r/L])  (✅ keep)

[
\boxed{,W^{(n)}(K);\propto;\left(\frac{L}{n}\right)^2
\left(1+\Big(\frac{K L}{n}\Big)^2\right)^{-3/2}
=\left(\frac{L}{n}\right)^2\left(1+\Big(\frac{K'}{n}\Big)^2\right)^{-3/2}.}
]
Your implementation matches this.

### 4.3 “1.5-power / transformed exponential”

If the intent is (\rho(r)=\exp[-(r/L)^{3/2}]), then the **n-fold similarity law** is
[
\rho(r)^n=\exp!\left[-\left(\frac{r}{L,n^{-2/3}}\right)^{3/2}\right]
\quad\Longrightarrow\quad
\boxed{,W^{(n)}(K)=L^2,n^{-4/3};\Phi!\big(K L;n^{-2/3}\big),}
]
for a **single** shape function (\Phi) (independent of (n)).
**Therefore:** in any analytic surrogate, the **special-function order must not depend on (n)**; only the **argument scales** ((n^{-2/3})) and the **amplitude** scales ((n^{-4/3})).

**Bug in MATLAB (must fix):** current code uses a modified-Bessel form with **order (1.5,n-1)**. That violates the similarity law and is not physically consistent.

**Two robust options:**

* **Numerical Hankel transform** (exact & stable):
  [
  W^{(n)}(K)=2\pi\int_0^\infty r,\exp!\Big[-n,(r/L)^{3/2}\Big],J_0(Kr),dr.
  ]
  Evaluate with Gauss–Laguerre (or mapped Gauss–Jacobi) for (r) and standard quadrature for (J_0).
* **Similarity-correct surrogate** (simple & fast):
  [
  \boxed{,W^{(n)}(K);\approx;\left(\frac{L}{n}\right)^2
  \left(1+\alpha^2\Big(\frac{K L}{n^{2/3}}\Big)^2\right)^{-3/2}},
  ]
  with (\alpha) chosen to match the curvature at (K=0). This preserves the required (n^{-4/3}) amplitude and (n^{-2/3}) argument scaling.

---

## 5) Single-scattering series (master formula and blocks)

Let (k_{iz}=k,cs), (k_{sz}=k,css).

### 5.1 Master formula

[
\boxed{,\sigma^0_{qp}
=\frac{k^2}{2};e^{-\sigma^2\left(k_{iz}^2+k_{sz}^2\right)}
\sum_{n=1}^{\infty}\frac{\sigma^{2n}}{n!};\big|I^{(n)}_{qp}\big|^2;W^{(n)}(K),,
\qquad q\in{h,v},;p\in{h,v}.}
]

### 5.2 The (I^{(n)}_{qp}) building block (exact structure)

[
\boxed{,I^{(n)}*{qp}
=(k*{sz}+k_{iz})^n,f_{qp},e^{-\sigma^2 k_{iz}k_{sz}}
+\frac{1}{4}\sum_{\beta=1}^{8};\mathcal{C}*\beta^{(n)};\mathcal{F}*{\beta};e^{-\sigma^2\Phi_\beta}, .}
]

* (f_{qp}): Kirchhoff field coefficients (your `fvv,fhh,fhv,fvh`), built from geometry, slopes (z_x,z_y), and the **transitioned** Fresnels (R_p^{(T)}) from §3.
* (\mathcal{F}_{\beta}): the eight complementary coefficients (your `fa..` for upper medium and `fb..` for lower medium).
* (\mathcal{C}_\beta^{(n)}): the power factors ((css\pm cs)^n,\ (css\pm q_t)^n,\ (cs\pm q_t)^n) appearing exactly as in your current `Ivv/Ihh/Ihv/Ivh` sums.
* (e^{-\sigma^2\Phi_\beta}): the Gaussian exponents. Your helper
  [
  \exp!\Big(-ks^2,[q^2 - q,(css - cs)]\Big)
  ]
  is fine **if and only if** you pass the right (q\in{cs,\ css,\ q_t}) for each of the eight terms. Keep that mapping one-to-one.

**Note on your long `fa*/fb*` formulas:** they match the required algebraic pattern (dependence on (z_x,z_y), ((1\pm R_{p})) factors, etc.). No structural change is required there.

### 5.3 Kirchhoff-only resummation (your “k-term”)

The Kirchhoff piece alone yields
[
\Big|(k_{sz}+k_{iz})^n f_{qp},e^{-\sigma^2 k_{iz}k_{sz}}\Big|^2
=|f_{qp}|^2;(k_{sz}+k_{iz})^{2n};e^{-2\sigma^2 k_{iz}k_{sz}}.
]
Your re-summation
[
\frac{1}{2},|f_{qp}|^2;\exp!\big(-ks^2(css+cs)^2\big);\sum_{n\ge 1}\left(\frac{ks^2(css+cs)^2}{n}\right)^{!n} W^{(n)}(K)
]
is consistent with the master formula (the loop’s running product by (1/n) implements (1/n!)).

---

## 6) Polarimetric reciprocity (monostatic symmetry)

For isotropic scalar media in **monostatic** geometry ((\theta_s=\theta_i), (\phi_s=\phi_i+\pi)):
[
\boxed{;\sigma^0_{hv}=\sigma^0_{vh};.}
]
Your code computes both; after the fixes above, they should agree numerically. For robustness, average them in monostatic use and warn if their relative difference exceeds a small tolerance (e.g., (10^{-6})).

---

## 7) Concrete defect list (as found) and exact fixes

1. **Specular half-angle “bug” (status update).**
   **Final verdict:** **no bug**. Your formula
   (csl=\sqrt{(1 + cs\cdot css - si\cdot sis\cdot csfs)/2})
   is correct because the half-angle is between **(-\hat{\mathbf k}_i)** and (\hat{\mathbf k}_s).

2. **Fresnel branch for lossy media (fix required).**
   Use (k_{tz}/k=\sqrt{\varepsilon_r-\sin^2\theta}) with the **decaying branch**: if (\Im{k_{tz}}<0) flip the sign. Then compute (R_h,R_v) from §2.2.

3. **Normal-incidence constants (fix required).**
   Set (r_0=(1-\sqrt{\varepsilon_r})/(1+\sqrt{\varepsilon_r})) and use **the same (r_0)** for both pols whenever a “normal-incidence” constant is used. The current code sets opposite signs; correct that.

4. **Transition function (replace legacy).**
   Implement (R_p^{(T)}=R_p(\theta_i)+[R_p(\theta_{\mathrm{sp}})-R_p(\theta_i)](1-S_p/S_p^{%280%29})).
   If any legacy block is temporarily kept:

   * Fix the **pol typo**: in the H-path use `rh0` not `rv0`.
   * Make Gaussian factors **consistent**: multiply by (\exp{-(ks\cdot cs)^2}) in **all** places; never divide by (\exp{+(ks\cdot cs)^2}).

5. **“1.5-power” spectrum (fix required).**
   Remove the Bessel-(K) form with order (1.5,n-1). Replace with either the **numerical Hankel** integral or the **similarity-correct surrogate** in §4.3.

6. **Complex near-singularity guards (small hygiene fix).**
   Replace tests like `abs(real(css-qslp))<tol` by `abs(css-qslp)<tol` (complex magnitude).

7. **Bessel symmetry (if any remains).**
   Prefer (K_{|\nu|}(x)) numerically (`besselk(abs(nu),x)`), since (K_{-\nu}=K_{\nu}).

---

## 8) Python port checklist (minimal, deterministic)

1. Inputs: (\theta_i,\theta_s,\phi_s,\ ks,\ kl,\ \varepsilon_r,) surface type (\in){Gaussian, Exponential, 1.5-power}.
2. Build shorthands (si,cs,sis,css,csfs,sfs).
3. Compute (K) and (csl=\cos\theta_{\mathrm{sp}}) using §0–§1.
4. Fresnels: compute (k_{tz}/k) with branch rule; evaluate (R_h(\theta_i),R_v(\theta_i)) and (R_h(\theta_{\mathrm{sp}}),R_v(\theta_{\mathrm{sp}})).
5. Transition: compute (R_p^{(T)}) by the (S_p/S_p^{(0)}) method.
6. Field coefficients: evaluate `fvv,fhh,fhv,fvh` and the eight complementary `fa*/fb*` with (R_p^{(T)}).
7. Spectra: (W^{(n)}(K)) per §4 (with corrected 1.5-power).
8. Series: assemble (I^{(n)}*{qp}), then (\sigma^0*{qp}) via §5.
9. (Monostatic) enforce/average HV and VH.
10. Return (\sigma^0) in linear or dB as needed.

---

## 9) Quick validation after patch

* **Small roughness** (ks\to 0): (\gamma_p) finite; (S_p^{(0)}) stable when computed with (ks=10^{-7}) vs (10^{-6}).
* **Angles**: monostatic ((\theta_s=\theta_i,\ \phi_s=\phi_i+\pi)) gives (\sigma^0_{hv}\approx\sigma^0_{vh}) and (\theta_{\mathrm{sp}}=0^\circ).
* **Spectra scaling**:
  Gaussian—width (\propto \sqrt{n}); Exponential—amplitude (\propto 1/n^2);
  1.5-power—curves collapse when plotting (n^{4/3}W^{(n)}) vs (K,L,n^{-2/3}).
* **Lossy soil**: (|R_p(\theta)|\le 1), smooth across near-critical/grazing angles (branch fixed).

---

### Final word

* Your **core AIEM machinery** (geometry, (K), (I^{(n)}) layout, Kirchhoff resummation, Gaussian/Exponential spectra) is correct.
* Apply these **decisive fixes**:

  * Fresnel **branch** & identical normal-incidence constant for both pols,
  * Transition function via (\gamma_p=1-S_p/S_p^{(0)}),
  * Replace the **1.5-power** spectrum with a similarity-correct method.

With these in place, the MATLAB implementation is physically sound and ready for an accurate, clean Python port.
