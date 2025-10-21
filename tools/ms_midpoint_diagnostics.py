#!/usr/bin/env python
"""Diagnostic script to inspect AIEM multiple-scattering midpoints.

This utility evaluates the “middle” LUT configuration (median RMS height) for
each correlation-length ratio in the 40° NMM3D table.  It records:

* Fraction of negative cells in the total multiple scattering integrand
* Kirchhoff (σ_kc) and complementary (σ_c) contributions
* Individual Kirchhoff term magnitudes (K1, K2, K3)
* Leading complementary term magnitudes (gc1–gc4) for quick triage

The goal is to make it easy to spot when destructive interference is being
suppressed (σ_kc overwhelmingly positive) or when complementary terms are
driving the integrand negative.

Usage
-----
    PYTHONPATH=src MPLCONFIGDIR=/tmp python tools/ms_midpoint_diagnostics.py
"""

from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

DATA_PATH = Path("data/NMM3D_LUT_NRCS_40degree.dat")
INCIDENCE_DEG = 40.0
NP_POINTS = 65
NMAX = 20


def _load_midpoint_rows() -> Iterable[Tuple[np.ndarray, float]]:
    """Yield the midpoint LUT row for each ℓ/σ ratio at 40° incidence."""
    lut = np.loadtxt(DATA_PATH)
    rows = lut[np.isclose(lut[:, 0], INCIDENCE_DEG)]
    ratios = sorted(set(rows[:, 1]))
    for ratio in ratios:
        subset = rows[np.isclose(rows[:, 1], ratio)]
        subset = subset[np.argsort(subset[:, 4])]  # sort by rms_norm
        yield subset[len(subset) // 2], ratio


def _integrate(arr: np.ndarray, weights: np.ndarray, real_only: bool = True) -> complex:
    """Weighted integral over the quadrature grid."""
    if real_only:
        arr = np.real(arr)
    return np.sum(arr * weights)


def main() -> None:
    # Local import so this tool does not run on machines without the module.
    from mwrtms.scattering.surface.iem import multiple_scattering as ms
    from mwrtms.scattering.surface.iem.fresnel_utils import compute_fresnel_incident
    from mwrtms.scattering.surface.iem.kirchhoff import compute_kirchhoff_coefficients

    diag = OrderedDict()
    for row, ratio in _load_midpoint_rows():
        (
            _theta_deg,
            _ratio,
            eps_r_real,
            eps_r_imag,
            rms_norm,
            _vv_db,
            _hh_db,
            _hv_db,
        ) = row

        frequency_ghz = 5.405
        lam = 0.3 / frequency_ghz
        k = 2 * np.pi / lam
        sigma = rms_norm * lam
        corr_length = ratio * sigma
        ks = k * sigma
        kl = k * corr_length

        geom = ms._prepare_geometry_params(
            np.deg2rad(INCIDENCE_DEG),
            np.deg2rad(INCIDENCE_DEG),
            0.0,
            np.deg2rad(180.0),
            k,
        )
        phys = ms.PhysicsParams(k=k, er=complex(eps_r_real, eps_r_imag))
        surf = ms.SurfaceParams("exponential", ks=ks, kl=kl, sigma=sigma)
        quad = ms._build_quadrature(surf, n_points=NP_POINTS, nmax=NMAX)

        U, V = quad.U, quad.V
        q1, q2 = ms._compute_vertical_wavenumbers(U, V, k, phys.er)
        weights = np.outer(quad.wu, quad.wv)
        mask = np.ones_like(U, dtype=bool)
        constants = ms._precompute_constants(surf, quad.Nmax)
        wn_provider = ms._make_Wn_provider(surf, constants)

        Rv, Rh, _ = compute_fresnel_incident(phys.er, geom.theta_i)
        fvv, fhh, fhv, fvh = compute_kirchhoff_coefficients(
            Rv,
            Rh,
            k,
            geom.theta_i,
            geom.theta_s,
            geom.phi_s,
            geom.phi_i,
        )
        f_map = {"vv": fvv, "hh": fhh, "hv": fhv, "vh": fvh}

        entry: Dict[str, Dict[str, float]] = {}
        propagator_diag: Dict[str, Dict[str, float]] = {}
        kernel_diag: Dict[str, Dict[str, float]] = {}
        for pol in ("vv", "hh", "hv"):
            propagators = ms._build_propagators(U, V, q1, q2, k, phys.er, geom, pol, False)
            propagator_diag[pol] = {
                "Fp_real_mean": float(np.real(propagators["Fp"]).mean()),
                "Fp_real_min": float(np.real(propagators["Fp"]).min()),
                "Fp_real_max": float(np.real(propagators["Fp"]).max()),
                "Fm_real_mean": float(np.real(propagators["Fm"]).mean()),
                "Gp_real_mean": float(np.real(propagators["Gp"]).mean()),
            }

            integrand_kc, integrand_c = ms._assemble_integrands_corrected(
                U,
                V,
                q1,
                q2,
                k,
                phys.er,
                geom,
                surf,
                wn_provider,
                quad.Nmax,
                pol,
                f_map.get(pol, 1.0 + 0.0j),
                constants,
                False,
                None,
            )

            total_real = np.real(integrand_kc + integrand_c)
            neg_frac = float(np.mean(total_real < 0))

            K1 = ms._build_gkc1(U, V, geom, q1, surf, wn_provider, quad.Nmax, constants, False)
            K2 = ms._build_gkc2(U, V, geom, q1, surf, wn_provider, quad.Nmax, constants, False)
            K3 = ms._build_gkc3(U, V, geom, q1, surf, wn_provider, quad.Nmax, constants, False)
            kernel_diag[pol] = {
                "K1_mean": float(np.real(K1).mean()),
                "K2_mean": float(np.real(K2).mean()),
                "K3_mean": float(np.real(K3).mean()),
            }

            kc_components = [
                np.conjugate(f_map[pol]) * _integrate(propagators["Fp"] * K1, weights, real_only=False),
                np.conjugate(f_map[pol]) * _integrate(propagators["Fm"] * K2, weights, real_only=False),
                np.conjugate(f_map[pol]) * _integrate(propagators["Gp"] * K3, weights, real_only=False),
            ]
            kc_sum = sum(kc_components)

            def shifted(name: str, idx: int | None) -> np.ndarray:
                return ms._evaluate_propagator_at_shifted_point(
                    propagators,
                    name,
                    idx,
                    U,
                    V,
                    k,
                    phys.er,
                    geom,
                    pol,
                    enable_guardrails=False,
                )

            C1 = ms._build_gc_block1(U, V, geom, q1, q1, surf, wn_provider, quad.Nmax, constants, False)
            C2 = ms._build_gc_block2(U, V, geom, q1, q1, surf, wn_provider, quad.Nmax, constants, False)
            comp_terms = [
                _integrate(propagators["Fp"] * np.conjugate(shifted("Fp", None)) * C1["gc1"], weights, real_only=False),
                _integrate(propagators["Fp"] * np.conjugate(shifted("Fm", 2)) * C1["gc2"], weights, real_only=False),
                _integrate(propagators["Fm"] * np.conjugate(shifted("Fp", 3)) * C1["gc3"], weights, real_only=False),
                _integrate(propagators["Fm"] * np.conjugate(shifted("Fm", 4)) * C1["gc4"], weights, real_only=False),
                _integrate(propagators["Gp"] * np.conjugate(shifted("Gp", 5)) * C1["gc5"], weights, real_only=False),
                _integrate(propagators["Gp"] * np.conjugate(shifted("Gm", 6)) * C1["gc6"], weights, real_only=False),
                _integrate(propagators["Gm"] * np.conjugate(shifted("Gp", 7)) * C1["gc7"], weights, real_only=False),
                _integrate(propagators["Gm"] * np.conjugate(shifted("Gm", 8)) * C1["gc8"], weights, real_only=False),
                _integrate(propagators["Fp"] * np.conjugate(shifted("Gp", 9)) * C2["gc9"], weights, real_only=False),
                _integrate(propagators["Fp"] * np.conjugate(shifted("Gm", 10)) * C2["gc10"], weights, real_only=False),
                _integrate(propagators["Fm"] * np.conjugate(shifted("Gp", 11)) * C2["gc11"], weights, real_only=False),
                _integrate(propagators["Fm"] * np.conjugate(shifted("Gm", 12)) * C2["gc12"], weights, real_only=False),
                _integrate(propagators["Gp"] * np.conjugate(shifted("Fp", 13)) * C2["gc13"], weights, real_only=False),
                _integrate(propagators["Gm"] * np.conjugate(shifted("Fp", 14)) * C2["gc14"], weights, real_only=False),
            ]
            comp_sum = sum(comp_terms)

            def serialize_terms(terms):
                return [
                    {
                        "real": float(np.real(term)),
                        "imag": float(np.imag(term)),
                    }
                    for term in terms
                ]

            entry[pol] = {
                "neg_frac": neg_frac,
                "sigma_kc": (k ** 2 / (8 * np.pi)) * float(np.real(kc_sum)),
                "sigma_comp": (k ** 2 / (64 * np.pi)) * float(np.real(comp_sum)),
                "sigma_ms": (k ** 2 / (8 * np.pi)) * float(np.real(kc_sum))
                + (k ** 2 / (64 * np.pi)) * float(np.real(comp_sum)),
                "kc_terms": serialize_terms(kc_components),
                "comp_terms_first4": serialize_terms(comp_terms[:4]),
            }

        diag[f"{ratio:g}"] = {
            pol.upper(): {
                "neg_frac": stats["neg_frac"],
                "sigma_kc": stats["sigma_kc"],
                "sigma_comp": stats["sigma_comp"],
                "sigma_ms": stats["sigma_ms"],
                "kc_terms": stats["kc_terms"],
                "comp_terms_first4": stats["comp_terms_first4"],
                "propagators": propagator_diag[pol],
                "kernels": kernel_diag[pol],
            }
            for pol, stats in entry.items()
        }

    print(json.dumps(diag, indent=2, default=float))


if __name__ == "__main__":
    main()
