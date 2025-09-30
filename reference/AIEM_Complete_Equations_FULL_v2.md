# AIEM — Complete Equations (Main Text 1–17) + Full Appendices (A1–A17, B1–B6, C1–C6 & B1–B6)

_Formatting_: Markdown + LaTeX math (for use in VS Code / Codex). Variables follow the paper’s notation.
- \(k\): free‑space wavenumber; \(\varepsilon_r, \mu_r\): relative permittivity, permeability.
- \(\mathbf{k}_i=(k_x,k_y,-q_1)\), \(\mathbf{k}_s=(k_{sx},k_{sy},k_{sz})\), with \(q_1=\sqrt{k^2-u^2-v^2}\), \(q_2=\sqrt{\varepsilon_r k^2-u^2-v^2}\).
- \(W^{(n)}(\cdot)\): roughness spectrum associated with the \(n\)-th power of the correlation function.
- \(R_v, R_h\): Fresnel reflection (upper interface); \(R=\tfrac{1}{2}(R_v-R_h)\) where used in cross‑pol blocks.

---

## Main equations

### (1) Decomposition of scattered field
$$
E_{qp}^{\,s}=E_{qp}^{\,k}+E_{qp}^{\,c}\,.
\tag{1}
$$

### (2) Kirchhoff (single‑bounce) field
$$
E_{qp}^{\,k}=C\,E_0\,f_{qp}\,e^{\,j(\mathbf{k}_s-\mathbf{k}_i)\cdot\mathbf{r}}
\quad\text{with}\quad
C=-\,\frac{jk}{4\pi R}e^{-jkR}.
\tag{2}
$$

### (3) Vectors
$$
\mathbf{k}_i=(k_x,k_y,-q_1),\qquad
\mathbf{k}_s=(k_{sx},k_{sy},k_{sz}).
\tag{3}
$$

### (4) Complementary field integral
$$
E^{\,c}_{qp}
=\frac{C\,E_0}{8\pi^2}\!\!\iiint\!\!\int
\!\left[
F_{qp}^{\pm}\,e^{j(u(x-x')+v(y-y'))-jq_1|z-z'|}
+
G_{qp}^{\pm}\,e^{j(u(x-x')+v(y-y'))-jq_2|z-z'|}
\right]
e^{\,j(\mathbf{k}_s\cdot\mathbf{r}-\mathbf{k}_i\cdot\mathbf{r}')} \,
dx\,dy\,dx'\,dy'\,du\,dv.
\tag{4}
$$

### (5) Upward/downward designation
$$
\{F_{qp},G_{qp}\}=
\begin{cases}
\{F^+_{qp},G^+_{qp}\}, & z>z',\\
\{F^-_{qp},G^-_{qp}\}, & z<z'.
\end{cases}
\tag{5}
$$

### (6) Mean scattered power
$$
P_s=\big\langle E^{\,s}_{qp}-\langle E^{\,s}_{qp}\rangle\big\rangle
\big\langle E^{\,s*}_{qp}-\langle E^{\,s*}_{qp}\rangle\big\rangle.
\tag{6}
$$

### (7) NRCS split
$$
\sigma^0_{qp}=\sigma_{qp}^{(s)}+\sigma_{qp}^{(m)}.
\tag{7}
$$

### (8) Single‑scattering split
$$
\sigma_{qp}^{(s)}=\sigma^{k}_{qp}+\sigma^{kc}_{qp}+\sigma^{c}_{qp}.
\tag{8}
$$

### (9) Explicit single‑scattering form
$$
\sigma_{qp}^{(s)}
=\frac{k^2}{2}\,e^{-\sigma^2(k_{sz}^2+k_z^2)}
\sum_{n=1}^{\infty}\frac{\sigma^{2n}}{n!}\,
\left|I^{n}_{qp}\right|^2\,
W^{(n)}(k_{sx}-k_x,\,k_{sy}-k_y),
\tag{9}
$$

with
$$
\begin{aligned}
I^n_{qp}=(k_{sz}+k_z)^n f_{qp}\,e^{-\sigma^2k_z k_{sz}}
+\frac{1}{4}\Big\{&
F^+_{qp}(-k_x,-k_y)\,(k_{sz}-k_z)^n e^{-\sigma^2(k_z^2-k_z k_{sz}+k_z^2)} \\
&+G^+_{qp}(-k_x,-k_y)\,(k_{sz}-k_{tz})^n e^{-\sigma^2(k_{tz}^2-k_{tz}k_{sz}+k_{tz}k_z)} \\
&+F^-_{qp}(-k_x,-k_y)\,(k_{sz}+k_z)^n e^{-\sigma^2(k_z^2+k_z k_{sz}-k_z^2)} \\
&+G^-_{qp}(-k_x,-k_y)\,(k_{sz}+k_{tz})^n e^{-\sigma^2(k_{tz}^2+k_{tz}k_{sz}-k_{tz}k_z)} \\
&+F^+_{qp}(-k_{sx},-k_{sy})\,(k_z+k_{sz})^n e^{-\sigma^2(k_{sz}^2-k_{sz}k_{sz}+k_{sz}k_z)} \\
&+G^+_{qp}(-k_{sx},-k_{sy})\,(k_z+k_{tsz})^n e^{-\sigma^2(k_{tsz}^2-k_{tsz}k_{sz}+k_{tsz}k_z)} \\
&+F^-_{qp}(-k_{sx},-k_{sy})\,(k_z-k_{sz})^n e^{-\sigma^2(k_{sz}^2+k_{sz}k_{sz}-k_{sz}k_z)} \\
&+G^-_{qp}(-k_{sx},-k_{sy})\,(k_z-k_{tsz})^n e^{-\sigma^2(k_{tsz}^2+k_{tsz}k_{sz}-k_{tsz}k_z)}
\Big\}.
\end{aligned}
\tag{10}
$$

### (11) Multiple‑scattering split (up to double‑bounce)
$$
\sigma_{qp}^{(m)}
=\sum_{l=1}^{3}\sigma^{kc_l}_{qp}(m)
+\sum_{i=1}^{8}\sigma^{c_i}_{qp}(m)
+\sum_{j=9}^{14}\sigma^{c_j}_{qp}(m).
\tag{11}
$$

### (12) Cross term (Kirchhoff × complementary)
$$
\sigma^{kc_l}_{qp}(m)
=\frac{k^2}{8\pi}\,\Re\!\left[
f_{qp}^* \iint\!\big\{
F^+_{qp}(u,v)\,g^{kc_l}(u,v,q_1)
+F^-_{qp}(u,v)\,g^{kc_l}(u,v,-q_1) \right.\\
\left.
+\,G^+_{qp}(u,v)\,g^{kc_l}(u,v,q_2)
+G^-_{qp}(u,v)\,g^{kc_l}(u,v,-q_2)
\big\}\,du\,dv\right].
\tag{12}
$$

### (13a) Complementary–complementary (part I)
$$
\begin{aligned}
\sigma^{c_i}_{qp}(m)=\frac{k^2}{64\pi}\iint\Big\{&
F^+_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_i(u,v;q_1,q'_1)
+F^+_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_i(u,v;q_1,-q'_1)\\
&+F^-_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_i(u,v;-q_1,q'_1)
+F^-_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_i(u,v;-q_1,-q'_1)\\
&+F^+_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_i(u,v;q_1,q'_2)
+F^+_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_i(u,v;q_1,-q'_2)\\
&+F^-_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_i(u,v;-q_1,q'_2)
+F^-_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_i(u,v;-q_1,-q'_2)\\
&+G^+_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_i(u,v;q_2,q'_1)
+G^+_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_i(u,v;q_2,-q'_1)\\
&+G^-_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_i(u,v;-q_2,q'_1)
+G^-_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_i(u,v;-q_2,-q'_1)\\
&+G^+_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_i(u,v;q_2,q'_2)
+G^+_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_i(u,v;q_2,-q'_2)\\
&+G^-_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_i(u,v;-q_2,q'_2)
+G^-_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_i(u,v;-q_2,-q'_2)
\Big\}\,du'\,dv'.
\end{aligned}
\tag{13a}
$$

### (13b) Complementary–complementary (part II, swapped arguments)
$$
\begin{aligned}
\sigma^{c_j}_{qp}(m)=\frac{k^2}{64\pi}\iint\Big\{&
F^+_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_j(u',v';q_1,q'_1)
+F^+_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_j(u',v';q_1,-q'_1)\\
&+F^-_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_j(u',v';-q_1,q'_1)
+F^-_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_j(u',v';-q_1,-q'_1)\\
&+F^+_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_j(u',v';q_1,q'_2)
+F^+_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_j(u',v';q_1,-q'_2)\\
&+F^-_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_j(u',v';-q_1,q'_2)
+F^-_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_j(u',v';-q_1,-q'_2)\\
&+G^+_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_j(u',v';q_2,q'_1)
+G^+_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_j(u',v';q_2,-q'_1)\\
&+G^-_{qp}(u,v)\,F^{+*}_{qp}(u',v')\,g^c_j(u',v';-q_2,q'_1)
+G^-_{qp}(u,v)\,F^{-*}_{qp}(u',v')\,g^c_j(u',v';-q_2,-q'_1)\\
&+G^+_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_j(u',v';q_2,q'_2)
+G^+_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_j(u',v';q_2,-q'_2)\\
&+G^-_{qp}(u,v)\,G^{+*}_{qp}(u',v')\,g^c_j(u',v';-q_2,q'_2)
+G^-_{qp}(u,v)\,G^{-*}_{qp}(u',v')\,g^c_j(u',v';-q_2,-q'_2)
\Big\}\,du'\,dv'.
\end{aligned}
\tag{13b}
$$

### (14) Argument mapping for \(g^c_i, g^c_j\)
$$
\begin{aligned}
&u'=u,\; v'=v,\; q'=q; && i=1,\\
&u'=-k_x-k_{sx}-u,\; v'=-k_y-k_{sy}-v; && i=2,\\
&u'=-k_{sx},\; v'=-k_{sy},\; q'=k_{sz}; && i=3,4,5,\\
&u'=-k_x,\; v'=-k_y,\; q'=k_z; && i=6,7,8,\\[2mm]
&u=-k_{sx},\; v=-k_{sy},\; q=k_{sz}; && j=9,10,11,\\
&u=-k_x,\; v=-k_y,\; q=k_z; && j=12,13,14,\\[2mm]
&q_1=\sqrt{k^2-u^2-v^2},\qquad & q'_1=\sqrt{k^2-u'^2-v'^2},\\
&q_2=\sqrt{\varepsilon_r k^2-u^2-v^2},\qquad & q'_2=\sqrt{\varepsilon_r k^2-u'^2-v'^2}.
\end{aligned}
\tag{14}
$$

### (15a) Slightly rough surface limit \((k\sigma\to 0)\)
$$
\sigma^{(m)}_{vh}=\sigma^{(m)}_{hv}
=\frac{2k^4\sigma^4\cos^2\theta}{\pi}\!\iint\!
\chi\,u^2v^2\,W(u-k\sin\theta, v)\,W(u+k\sin\theta, v)\,du\,dv.
\tag{15a}
$$

### (15b) Definitions for (15a)
$$
\begin{aligned}
\chi &= \left|\gamma\right|^2 + \frac{d_0^*}{d_0}\,\gamma\kappa^* + \frac{d_0}{d_0^*}\,\gamma^*\kappa + \left|\frac{\kappa}{d_0}\right|^2,\\
\gamma &= \frac{2R^2}{q_1},\qquad R=R_v-R_h,\\
\kappa &= \frac{(\varepsilon_r^2+6\varepsilon_r+1)R^2 - 2(\varepsilon_r^2-1)R + (\varepsilon_r-1)^2}{4q_2\,\varepsilon_r},\\
d_0 &= q_2^2 - k^2\cos^2\theta,\qquad
q_1=\sqrt{k^2-u^2-v^2},\quad q_2=\sqrt{\varepsilon_r k^2-u^2-v^2}.
\end{aligned}
\tag{15b}
$$

### (16) SPM cross‑pol (for comparison)
$$
\sigma^{r}_{hv}=\sigma^{r}_{vh}=
\frac{2k^4\sigma^4\cos^2\theta}{\pi}\,
\frac{|(\varepsilon_r-1)(R_v-R_h)|^2}{|k'_z+\varepsilon_r k_z|^2}\,
\iint u^2v^2\,W(u-k\sin\theta, v)\,W(u+k\sin\theta, v)\,du\,dv,
\tag{16}
$$
with \(k_z=\sqrt{k^2-u^2-v^2}\), \(k'_z=\sqrt{\varepsilon_r k^2-u'^2-v'^2}\).

### (17) PEC limit
$$
\lim_{\varepsilon_r\to\infty}\sigma^{(m)}_{vh}
=\lim_{\varepsilon_r\to\infty}\sigma^{(m)}_{hv}
=\frac{8k^4\sigma^4\cos^2\theta}{\pi}\!\iint\!
\frac{u^2v^2}{q_1^2}\,W(u-k\sin\theta, v)\,W(u+k\sin\theta, v)\,du\,dv.
\tag{17}
$$

---

## Appendix A — Factors \(g_{kc\ell},\,g_i^c,\,g_j^c\)

Below, all exponentials share the common factor
\(e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}\)
with appropriate \(q\)–\(q'\) choices per equation; arguments of \(W^{(n)}\) are given explicitly.

### \(g_{kc\ell}(u,v,q)\)

$$
\begin{aligned}
g_{kc1}(u,v,q)
&=e^{-\sigma^2(k_{sz}^2+k_z^2+k_{sz}k_z+q^2-qk_{sz}+qk_z)}\,
\sum_{m=1}^{\infty}\frac{\left[\sigma^2(k_z+q)(k_{sz}+k_z)\right]^m}{m!}\,
W^{(m)}(k_x+u,\,k_y+v)\\
&\quad\times
\sum_{n=1}^{\infty}\frac{\left[\sigma^2(k_{sz}-q)(k_{sz}+k_z)\right]^n}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v).
\end{aligned}
\tag{A1}
$$

$$
\begin{aligned}
g_{kc2}(u,v,q)
&=e^{-\sigma^2(k_{sz}^2+k_z^2+k_{sz}k_z+q^2-qk_{sz}+qk_z)}\,
\sum_{m=1}^{\infty}\frac{\left[\sigma^2(k_z+q)(k_{sz}+k_z)\right]^m}{m!}\,
W^{(m)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times
\sum_{n=1}^{\infty}\frac{\left[-\sigma^2(k_{sz}-q)(k_z+q)\right]^n}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v).
\end{aligned}
\tag{A2}
$$

$$
\begin{aligned}
g_{kc3}(u,v,q)
&=e^{-\sigma^2(k_{sz}^2+k_z^2+k_{sz}k_z+q^2-qk_{sz}+qk_z)}\,
\sum_{m=1}^{\infty}\frac{\left[-\sigma^2(k_{sz}-q)(k_z+q)\right]^m}{m!}\,
W^{(m)}(k_x+u,\,k_y+v)\\
&\quad\times
\sum_{n=1}^{\infty}\frac{\left[\sigma^2(k_{sz}-q)(k_{sz}+k_z)\right]^n}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy}).
\end{aligned}
\tag{A3}
$$

### \(g_i^c(u,v;q,q')\)  (i = 1…8)

$$
\begin{aligned}
g_1^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v)\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_z+q')]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A4}
$$

$$
\begin{aligned}
g_2^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_z+q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v)\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_{sz}-q')]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A5}
$$

$$
\begin{aligned}
g_3^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_z+q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v)\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_z+q')]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A6}
$$

$$
\begin{aligned}
g_4^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_z+q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[-\sigma^2(k_{sz}-q)(k_z+q)]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A7}
$$

$$
\begin{aligned}
g_5^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_z+q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[-\sigma^2(k_{sz}-q)(k_z+q)]^{m}}{m!}\,
W^{(m)}(k_{sx}+u,\,k_{sy}+v).
\end{aligned}
\tag{A8}
$$

$$
\begin{aligned}
g_6^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u,\,k_{sy}+v)\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_z+q')]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A9}
$$

$$
\begin{aligned}
g_7^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[-\sigma^2(k_{sz}-q)(k_z+q)]^{m}}{m!}\,
W^{(m)}(k_x+u,\,k_y+v).
\end{aligned}
\tag{A10}
$$

$$
\begin{aligned}
g_8^c(u,v;q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[-\sigma^2(k_{sz}-q)(k_z+q)]^{m}}{m!}\,
W^{(m)}(k_{sx}+u,\,k_{sy}+v).
\end{aligned}
\tag{A11}
$$

### \(g_j^c(u',v';q,q')\)  (j = 9…14)

$$
\begin{aligned}
g_9^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u',\,k_{sy}+v')\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_z+q')]^{m}}{m!}\,
W^{(m)}(k_x+u',\,k_y+v').
\end{aligned}
\tag{A12}
$$

$$
\begin{aligned}
g_{10}^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\,-\sigma^2(k_{sz}-q')(k_z+q')\,]^{m}}{m!}\,
W^{(m)}(k_x+u',\,k_y+v').
\end{aligned}
\tag{A13}
$$

$$
\begin{aligned}
g_{11}^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_z+q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\,-\sigma^2(k_{sz}-q')(k_z+q')\,]^{m}}{m!}\,
W^{(m)}(k_{sx}+u',\,k_{sy}+v').
\end{aligned}
\tag{A14}
$$

$$
\begin{aligned}
g_{12}^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_{sx}+u',\,k_{sy}+v')\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_z+q')]^{m}}{m!}\,
W^{(m)}(k_x+u',\,k_y+v').
\end{aligned}
\tag{A15}
$$

$$
\begin{aligned}
g_{13}^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_{sz}-q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\,-\sigma^2(k_{sz}-q')(k_z+q')\,]^{m}}{m!}\,
W^{(m)}(k_x+u',\,k_y+v').
\end{aligned}
\tag{A16}
$$

$$
\begin{aligned}
g_{14}^c(u',v';q,q')
&=e^{-\sigma^2(k_{sz}^2+k_z^2+q^2+q'^2-(k_{sz}-k_z)(q+q'))}
\sum_{n=1}^{\infty}\frac{[\sigma^2(k_z+q)(k_{sz}-q')]^{n}}{n!}\,
W^{(n)}(k_x-k_{sx},\,k_y-k_{sy})\\
&\quad\times\sum_{m=1}^{\infty}\frac{[\,-\sigma^2(k_{sz}-q')(k_z+q')\,]^{m}}{m!}\,
W^{(m)}(k_{sx}+u',\,k_{sy}+v').
\end{aligned}
\tag{A17}
$$

---

## Appendix B — Upward reradiation coefficients \(F^+_{qp}, G^+_{qp}\)

(Downward ones \(F^-_{qp},G^-_{qp}\) are analogous by sign flips and \(q_1\!\leftrightarrow\! -q_1\), \(q_2\!\leftrightarrow\! -q_2\).)

### VV
$$
\begin{aligned}
F^+_{vv}&= -\Big(\frac{1-R_v}{q_1}\Big)(1+R_v)C_1
+\Big(\frac{1-R_v}{q_1}\Big)(1-R_v)C_2
+\Big(\frac{1-R_v}{q_1}\Big)(1+R_v)C_3\\
&\quad+\Big(\frac{1+R_v}{q_1}\Big)(1-R_v)C_4
+\Big(\frac{1+R_v}{q_1}\Big)(1+R_v)C_5
+\Big(\frac{1+R_v}{q_1}\Big)(1-R_v)C_6,\\[1mm]
G^+_{vv}&=\Big(\frac{(1+R_v)\,u_r}{q_2}\Big)(1+R_v)C_{1t}
-\Big(\frac{1+R_v}{q_2}\Big)(1-R_v)C_{2t}
-\Big(\frac{1+R_v\,\varepsilon_r}{q_2}\Big)(1+R_v)C_{3t}\\
&\quad-\Big(\frac{(1-R_v)\,\varepsilon_r}{q_2}\Big)(1-R_v)C_{4t}
-\Big(\frac{1-R_v}{q_2}\Big)(1+R_v)C_{5t}
-\Big(\frac{(1-R_v)\,u_r}{q_2}\Big)(1-R_v)C_{6t}.
\end{aligned}
\tag{B1–B2}
$$

### HH
$$
\begin{aligned}
F^+_{hh}&=\Big(\frac{1-R_h}{q_1}\Big)(1+R_h)C_1
-\Big(\frac{1-R_h}{q_1}\Big)(1-R_h)C_2
-\Big(\frac{1-R_h}{q_1}\Big)(1+R_h)C_3\\
&\quad-\Big(\frac{1+R_h}{q_1}\Big)(1-R_h)C_4
-\Big(\frac{1+R_h}{q_1}\Big)(1+R_h)C_5
-\Big(\frac{1+R_h}{q_1}\Big)(1-R_h)C_6,\\[1mm]
G^+_{hh}&=-\Big(\frac{(1+R_h)\varepsilon_r}{q_2}\Big)(1+R_h)C_{1t}
+\Big(\frac{1+R_h}{q_2}\Big)(1-R_h)C_{2t}
+\Big(\frac{(1+R_h)u_r}{q_2}\Big)(1+R_h)C_{3t}\\
&\quad+\Big(\frac{(1-R_h)u_r}{q_2}\Big)(1-R_h)C_{4t}
+\Big(\frac{1-R_h}{q_2}\Big)(1+R_h)C_{5t}
+\Big(\frac{(1-R_h)\varepsilon_r}{q_2}\Big)(1-R_h)C_{6t}.
\end{aligned}
\tag{B3–B4}
$$

### HV (cross‑pol)
$$
\begin{aligned}
F^+_{hv}&=\Big(\frac{1-R}{q_1}\Big)(1+R)B_1
-\Big(\frac{1-R}{q_1}\Big)(1-R)B_2
-\Big(\frac{1-R}{q_1}\Big)(1+R)B_3\\
&\quad+\Big(\frac{1+R}{q_1}\Big)(1-R)B_4
+\Big(\frac{1+R}{q_1}\Big)(1+R)B_5
+\Big(\frac{1+R}{q_1}\Big)(1-R)B_6,\\[1mm]
G^+_{hv}&=-\Big(\frac{(1+R)\mu_r}{q_2}\Big)(1+R)B_{1t}
+\Big(\frac{1+R}{q_2}\Big)(1-R)B_{2t}
+\Big(\frac{1+R\,\varepsilon_r}{q_2}\Big)(1+R)B_{3t}\\
&\quad-\Big(\frac{(1-R)\varepsilon_r}{q_2}\Big)(1-R)B_{4t}
-\Big(\frac{1-R}{q_2}\Big)(1+R)B_{5t}
-\Big(\frac{(1-R)\mu_r}{q_2}\Big)(1-R)B_{6t}.
\end{aligned}
\tag{B5–B6}
$$

---

## Appendix C — Coefficients \(C_i\) and \(B_i\) (take \(\varphi_i=0\))

Define \(\theta_s\equiv\theta\) for backscatter; \((u,v)\) are spectral components; \((z_x,z_y)\) and \((z'_x,z'_y)\) are direction cosines of incident and scattering unit vectors.

### \(C_1\)–\(C_6\)
$$
\begin{aligned}
C_1&=-\cos\varphi_s(-\cos\varphi - z_x z'_x\cos\varphi - z_x z'_y\sin\varphi)\\
&\quad+\sin\varphi_s(\sin\varphi + z'_x z_y\cos\varphi + z_y z'_y\sin\varphi),\\[1mm]
C_2&=-\cos\varphi_s(-q_1\cos\theta\cos\varphi - u z_x\cos\theta - v z'_y\cos\theta\cos\varphi - q_1 z'_x\sin\theta \\
&\qquad - u z_x z'_x\sin\theta - v z_x z'_y \sin\theta - v z_x\cos\theta\sin\varphi + v z'_x\cos\theta\sin\varphi) \\
&\quad+\sin\varphi_s(u z_y\cos\theta\cos\varphi - u z'_y\cos\theta + u z'_x z_y \sin\theta + q_1 z'_y\sin\theta + v z_y z'_y \sin\theta \\
&\qquad + q_1\cos\theta\sin\varphi + u z'_x\cos\theta\sin\varphi + v z_y\cos\theta\sin\varphi),\\[1mm]
C_3&=\cos\varphi_s(u z'_x\cos\theta\cos\varphi - q_1 z_x z'_x\cos\theta\cos\varphi - u\sin\theta + q_1 z_x\sin\theta \\
&\qquad + u z'_y\cos\theta\sin\varphi - q_1 z_x z'_y\cos\theta\sin\varphi)\\
&\quad+\sin\varphi_s(v z'_x\cos\theta\cos\varphi - q_1 z'_x z_y\cos\theta\cos\varphi - v\sin\theta + q_1 z_y\sin\theta \\
&\qquad + v z'_y\cos\theta\sin\varphi - q_1 z_y z'_y\cos\theta\sin\varphi),\\[1mm]
C_4&=\sin\theta_s(-z_x\cos\theta\cos\varphi - z_x z'_x\sin\theta - z_y z'_y\sin\theta - z_y\cos\theta\sin\varphi)\\
&\quad-\cos\theta_s\cos\varphi_s(-\cos\theta\cos\varphi - z_y z'_y\cos\theta\cos\varphi - z'_x\sin\theta + z'_x z_y\cos\theta\sin\varphi)\\
&\quad-\cos\theta_s\sin\varphi_s(z_x z'_y\cos\theta\cos\varphi - z'_y\sin\theta - \cos\theta\sin\varphi - z_x z'_x\cos\theta\sin\varphi),\\[1mm]
C_5&=\sin\theta_s(q_1 z_x\cos\varphi + u z_x z'_x\cos\varphi + v z'_x z_y\cos\varphi + q_1 z_y\sin\varphi + u z_x z'_y\sin\varphi + v z_y z'_y\sin\varphi)\\
&\quad-\cos\theta_s\cos\varphi_s(q_1\cos\varphi + u z'_x\cos\varphi + v z_y\cos\varphi - u z_y\sin\varphi + u z'_y\sin\varphi)\\
&\quad-\cos\theta_s\sin\varphi_s(-v z_x\cos\varphi + v z'_x\cos\varphi + q_1\sin\varphi + u z_x\sin\varphi + v z'_y\sin\varphi),\\[1mm]
C_6&=\sin\theta_s(v z_x z'_y\cos\varphi - u z_y z'_y\cos\varphi - v z_x z'_x\sin\varphi + u z'_x z_y\sin\varphi)\\
&\quad+\cos\theta_s\cos\varphi_s(-v z'_y\cos\varphi + q_1 z_y z'_y\cos\varphi + v z'_x\sin\varphi - q_1 z'_x z_y\sin\varphi)\\
&\quad+\cos\theta_s\sin\varphi_s(u z'_y\cos\varphi - q_1 z_x z'_y\cos\varphi - u z'_x\sin\varphi + q_1 z_x z'_x\sin\varphi).
\end{aligned}
$$

### \(B_1\)–\(B_6\)
$$
\begin{aligned}
B_1&=\sin\theta_s(-z_y\cos\varphi+z_x\sin\varphi)
-\cos\theta_s\cos\varphi_s(z'_x z_y\cos\varphi+\sin\varphi+z_y z'_y\sin\varphi)\\
&\quad-\cos\theta_s\sin\varphi_s(-\cos\varphi - z_x z'_x\cos\varphi - z_x z'_y\sin\varphi),\\[1mm]
B_2&=\sin\theta_s(-q_1 z_y\cos\theta\cos\varphi - u z_x z'_y\cos\theta\cos\varphi - v z_y z'_y\cos\theta\cos\varphi - q_1 z'_x z_y \sin\theta \\
&\qquad + q_1 z_x z'_y \sin\theta + q_1 z_x\cos\theta\sin\varphi + u z_x z'_x\cos\theta\sin\varphi + v z'_x z_y\cos\theta\sin\varphi)\\
&\quad-\cos\theta_s\cos\varphi_s(u z_y\cos\theta\cos\varphi - u z'_y\cos\theta\cos\varphi + u z'_x z_y\sin\theta + q_1 z'_y\sin\theta + v z_y z'_y\sin\theta \\
&\qquad + q_1\cos\theta\sin\varphi + u z'_x\cos\theta\sin\varphi + v z_y\cos\theta\sin\varphi)\\
&\quad-\cos\theta_s\sin\varphi_s(-q_1\cos\theta\cos\varphi - u z_x\cos\theta\cos\varphi - v z'_y\cos\theta\cos\varphi - q_1 z'_x\sin\theta \\
&\qquad - u z_x z'_x\sin\theta - v z_x z'_y\sin\theta - v z_x\cos\theta\sin\varphi + v z'_x\cos\theta\sin\varphi),\\[1mm]
B_3&=\sin\theta_s(v z_x z'_x\cos\theta\cos\varphi - u z_y z'_x\cos\theta\cos\varphi - v z_x\sin\theta + u z_y\sin\theta \\
&\qquad + v z_x z'_y\cos\theta\sin\varphi - u z_y z'_y\cos\theta\sin\varphi)\\
&\quad+\cos\theta_s\cos\varphi_s(-v z'_x\cos\theta\cos\varphi + q_1 z'_x z_y\cos\theta\cos\varphi + v\sin\theta - q_1 z_y\sin\theta \\
&\qquad - v z'_y\cos\theta\sin\varphi + q_1 z_y z'_y\cos\theta\sin\varphi)\\
&\quad+\cos\theta_s\sin\varphi_s(u z'_x\cos\theta\cos\varphi - q_1 z_x z'_x\cos\theta\cos\varphi - u\sin\theta + q_1 z_x\sin\theta \\
&\qquad + u z'_y\cos\theta\sin\varphi - q_1 z_x z'_y\cos\theta\sin\varphi),\\[1mm]
B_4&=-\cos\varphi_s(z_x z'_y\cos\theta\cos\varphi - z'_y\sin\theta - \cos\theta\sin\varphi - z_x z'_x\cos\theta\sin\varphi)\\
&\quad+\sin\varphi_s(-\cos\theta\cos\varphi - z_y z'_y\cos\theta\cos\varphi - z'_x\sin\theta + z'_x z_y\cos\theta\sin\varphi),\\[1mm]
B_5&=-\cos\varphi_s(-v z_x\cos\varphi + v z'_x\cos\varphi + q_1\sin\varphi + u z_x\sin\varphi - v z'_y\sin\varphi)\\
&\quad+\sin\varphi_s(q_1\cos\varphi + u z'_x\cos\varphi + v z_y\cos\varphi - u z_y\sin\varphi + u z'_y\sin\varphi),\\[1mm]
B_6&=\cos\varphi_s(u z'_y\cos\varphi - q_1 z_x z'_y\cos\varphi - u z'_x\sin\varphi + q_1 z_x z'_x\sin\varphi)\\
&\quad+\sin\varphi_s(v z'_y\cos\varphi - q_1 z_y z'_y\cos\varphi - v z'_x\sin\varphi + q_1 z'_x z_y\sin\varphi).
\end{aligned}
$$

> Coefficients with subscript \(t\) (e.g., \(C_{it}, B_{it}\)) for the **lower medium** are obtained by replacing \(q_1\) with \(q_2\).

---

### Notes
- All equations were transcribed verbatim (structure and terms) from the uploaded PDF. If you plan to plug these into code, prefer **symbolic functions** for \(W^{(n)}\) and supply a concrete spectrum (e.g., Gaussian or exponential) at runtime.
- Minor typographic harmonization was applied (e.g., consistent subscripts \(k_{sz}, k_z\); explicit \(W^{(n)}\) arguments).

