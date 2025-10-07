function [sigma_vv, sigma_hh, sigma_hv, sigma_vh] = aiem_single_scattering(theta_i, theta_s, phi_s, kl, ks, eps_r, correlation_type)
% AIEM_SINGLE_SCATTERING Corrected AIEM single scattering implementation
%
% Computes backscatter coefficients using the Advanced Integral Equation Model
% with all bug fixes applied from the Python implementation.
%
% Inputs:
%   theta_i          - Incident angle (degrees)
%   theta_s          - Scattered angle (degrees)
%   phi_s            - Scattered azimuth angle (degrees)
%   kl               - Normalized correlation length (k * L)
%   ks               - Normalized RMS height (k * sigma)
%   eps_r            - Complex relative permittivity
%   correlation_type - Surface correlation function:
%                      1 = Gaussian
%                      2 = Exponential
%                      3 = 1.5-power
%
% Outputs:
%   sigma_vv, sigma_hh, sigma_hv, sigma_vh - Backscatter coefficients (linear)
%
% Bug fixes applied:
%   1. Fresnel branch correction for lossy media
%   2. Normal incidence: rh0 = rv0 (not -rv0)
%   3. Transition function: use rh0 in H-path
%   4. 1.5-power spectrum: similarity-correct formula
%   5. Complex magnitude checks in complementary coefficients
%
% Note: This uses the LEGACY transition function. For <1 dB RMSE vs NMM3D,
%       implement the new S_p/S_p^(0) method from the bug report.

    % Convert angles to radians
    theta_i = deg2rad(theta_i);
    theta_s = deg2rad(theta_s);
    phi_s = deg2rad(phi_s);
    phi_i = 0.0;  % Standard convention
    
    % Trigonometric quantities
    si = sin(theta_i);
    cs = cos(theta_i);
    sis = sin(theta_s);
    css = cos(theta_s);
    sfs = sin(phi_s);
    csfs = cos(phi_s);
    
    % Compute spatial frequency K
    K = kl * sqrt((sis * csfs - si)^2 + (sis * sfs)^2);
    
    % Determine number of spectral terms
    n_terms = determine_spectral_terms(ks, cs, css);
    
    % Compute roughness spectrum for all orders
    spectra = zeros(1, n_terms);
    for n = 1:n_terms
        spectra(n) = roughness_spectrum(kl, K, n, correlation_type);
    end
    
    % Compute Fresnel coefficients (with bug fixes)
    [Rvi, Rhi, Rvhi, Rvl, Rhl, Rvhl, rv0, rh0] = fresnel_coefficients(eps_r, theta_i, theta_s, phi_s);
    
    % Compute transition function (with bug fixes)
    [Tfv, Tfh] = transition_function(eps_r, theta_i, ks, cs, spectra, n_terms);
    
    % Transitioned reflection coefficients
    Rvtran = Rvi + (Rvl - Rvi) * Tfv;
    Rhtran = Rhi + (Rhl - Rhi) * Tfh;
    Rvhtran = (Rvtran - Rhtran) / 2.0;
    
    % Compute Kirchhoff field coefficients
    k = ks / (ks / kl * kl);  % Reconstruct k from ks and kl (not used in coefficients)
    [fvv, fhh, fhv, fvh] = kirchhoff_coefficients(Rvtran, Rhtran, 1.0, theta_i, theta_s, phi_s, phi_i);
    
    % Compute Kirchhoff term for each polarization
    kterm_vv = compute_kirchhoff_term(fvv, ks, cs, css, spectra, n_terms);
    kterm_hh = compute_kirchhoff_term(fhh, ks, cs, css, spectra, n_terms);
    kterm_hv = compute_kirchhoff_term(fhv, ks, cs, css, spectra, n_terms);
    kterm_vh = compute_kirchhoff_term(fvh, ks, cs, css, spectra, n_terms);
    
    % Compute complementary term for each polarization
    cterm_vv = compute_complementary_term(Rvi, Rvhi, eps_r, ks, si, sis, cs, css, sfs, csfs, fvv, spectra, n_terms, 'vv');
    cterm_hh = compute_complementary_term(Rhi, Rvhi, eps_r, ks, si, sis, cs, css, sfs, csfs, fhh, spectra, n_terms, 'hh');
    cterm_hv = compute_complementary_term(Rvhi, Rvhi, eps_r, ks, si, sis, cs, css, sfs, csfs, fhv, spectra, n_terms, 'hv');
    cterm_vh = compute_complementary_term(Rvhi, Rvhi, eps_r, ks, si, sis, cs, css, sfs, csfs, fvh, spectra, n_terms, 'vh');
    
    % Total single scattering coefficients
    sigma_vv = real(kterm_vv + cterm_vv);
    sigma_hh = real(kterm_hh + cterm_hh);
    sigma_hv = real(kterm_hv + cterm_hv);
    sigma_vh = real(kterm_vh + cterm_vh);
end


function n_terms = determine_spectral_terms(ks, cs, css)
    % Determine number of spectral terms for series convergence
    tolerance = 1e-16;
    iterm = 1;
    temp_old = 0.0;
    temp = ks^2 * (cs + css)^2;
    
    while abs(temp - temp_old) > tolerance && iterm < 1000
        temp_old = temp;
        iterm = iterm + 1;
        temp = temp_old * (ks^2 * (cs + css)^2) / iterm;
    end
    
    n_terms = iterm;
end


function kterm = compute_kirchhoff_term(f, ks, cs, css, spectra, n_terms)
    % Compute Kirchhoff (specular) scattering term
    sum_val = 0.0;
    temp = 1.0;
    ks2 = ks^2;
    
    for n = 1:n_terms
        temp = temp * (ks2 * (cs + css)^2) / n;
        sum_val = sum_val + temp * spectra(n);
    end
    
    expk = exp(-ks2 * (css + cs)^2) * sum_val;
    kterm = 0.5 * expk * abs(f)^2;
end


function cterm = compute_complementary_term(R, Rhv, eps_r, ks, si, sis, cs, css, sfs, csfs, f, spectra, n_terms, pol)
    % Compute complementary (multiple scattering) term
    
    % Wave vector components
    qq = cs;
    qqt = sqrt(eps_r - si^2);
    qqs = css;
    qqts = sqrt(eps_r - sis^2);
    
    % Compute scattering integrals
    I_pol = compute_scattering_integral(R, Rhv, eps_r, ks, si, sis, cs, css, sfs, csfs, ...
                                        qq, qqt, qqs, qqts, f, n_terms, pol);
    
    % Compute final complementary term
    sum_val = 0.0;
    temp = 1.0;
    ks2 = ks^2;
    cs2 = cs^2;
    css2 = css^2;
    
    for n = 1:n_terms
        temp = temp * (ks2 / n);
        sum_val = sum_val + temp * abs(I_pol(n))^2 * spectra(n);
    end
    
    cterm = 0.5 * exp(-ks2 * (cs2 + css2)) * sum_val;
end


function I_pol = compute_scattering_integral(R, Rhv, eps_r, ks, si, sis, cs, css, sfs, csfs, ...
                                             qq, qqt, qqs, qqts, f, n_terms, pol)
    % Compute scattering integral for specified polarization
    I_pol = zeros(1, n_terms);
    ks2 = ks^2;
    
    % Wave vector components for complementary field
    qq1 = qq;
    qq2 = qqs;
    qq5 = qqt;
    qq6 = qqts;
    
    % Exponential helper function
    expal = @(q) exp(-ks2 * (q^2 - q * (css - cs)));
    
    % Compute complementary field coefficients for all 8 branches
    % Air-side, incident up
    Faupi = complementary_coefficients(-si, 0.0, qq1, qq1, qq, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, false) * expal(qq1);
    % Air-side, incident down
    Fadni = complementary_coefficients(-si, 0.0, -qq1, -qq1, qq, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, false) * expal(-qq1);
    % Air-side, scattered up
    Faups = complementary_coefficients(-sis*csfs, -sis*sfs, qq2, qq2, qqs, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, false) * expal(qq2);
    % Air-side, scattered down
    Fadns = complementary_coefficients(-sis*csfs, -sis*sfs, -qq2, -qq2, qqs, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, false) * expal(-qq2);
    % Substrate-side, incident up
    Fbupi = complementary_coefficients(-si, 0.0, qqt, qq5, qqt, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, true) * expal(qq5);
    % Substrate-side, incident down
    Fbdni = complementary_coefficients(-si, 0.0, -qqt, -qq5, qqt, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, true) * expal(-qq5);
    % Substrate-side, scattered up
    Fbups = complementary_coefficients(-sis*csfs, -sis*sfs, qqts, qq6, qqts, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, true) * expal(qq6);
    % Substrate-side, scattered down
    Fbdns = complementary_coefficients(-sis*csfs, -sis*sfs, -qqts, -qq6, qqts, R, eps_r, ...
                                       si, sis, cs, css, sfs, csfs, pol, true) * expal(-qq6);
    
    % Compute scattering integrals for each order
    for n = 1:n_terms
        I_pol(n) = (cs + css)^n * f * exp(-ks2 * cs * css) + ...
                   0.25 * (Faupi * (css - qq1)^n + ...
                           Fadni * (css + qq1)^n + ...
                           Faups * (cs + qq2)^n + ...
                           Fadns * (cs - qq2)^n + ...
                           Fbupi * (css - qq5)^n + ...
                           Fbdni * (css + qq5)^n + ...
                           Fbups * (cs + qq6)^n + ...
                           Fbdns * (cs - qq6)^n);
    end
end
