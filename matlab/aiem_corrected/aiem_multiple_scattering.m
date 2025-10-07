function sigma_ms = aiem_multiple_scattering(theta_i, theta_s, phi_s, kl, ks, eps_r, sigma, correlation_type, polarization, varargin)
% AIEM_MULTIPLE_SCATTERING Complete multiple scattering implementation for AIEM
%
% This is a COMPLETE translation of the Python implementation without simplifications.
% Implements second-order multiple scattering following Yang et al. (2017).
%
% USAGE:
%   sigma_ms = aiem_multiple_scattering(theta_i, theta_s, phi_s, kl, ks, eps_r, ...
%                                       sigma, correlation_type, polarization)
%   sigma_ms = aiem_multiple_scattering(..., 'n_points', 129, 'nmax', 8)
%
% INPUTS:
%   theta_i          - Incident angle (radians)
%   theta_s          - Scattered angle (radians)
%   phi_s            - Scattered azimuth (radians)
%   kl               - Normalized correlation length (k * L)
%   ks               - Normalized RMS height (k * sigma)
%   eps_r            - Complex relative permittivity
%   sigma            - RMS height (meters)
%   correlation_type - Surface correlation: 1=Gaussian, 2=Exponential
%   polarization     - Polarization: 'vv', 'hh', 'hv', 'vh'
%
% OPTIONAL NAME-VALUE PAIRS:
%   'n_points' - Number of quadrature points per dimension (default: 129)
%   'nmax'     - Maximum order for spectral series (default: 8)
%
% OUTPUT:
%   sigma_ms - Multiple scattering coefficient (linear)
%
% REFERENCE:
%   Yang, Y., Chen, K. S., Tsang, L., & Yu, L. (2017). "Depolarized
%   backscattering of rough surface by AIEM model." IEEE Journal of Selected
%   Topics in Applied Earth Observations and Remote Sensing, 10(11), 4740-4752.
%
% AUTHOR: Complete translation from Python implementation
% DATE: 2024

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'n_points', 129, @isnumeric);
    addParameter(p, 'nmax', 8, @isnumeric);
    parse(p, varargin{:});
    n_points = p.Results.n_points;
    nmax = p.Results.nmax;
    
    % Wavenumber
    k = ks / sigma;
    phi_i = 0.0;  % Standard convention
    
    % Prepare geometry parameters
    geom = prepare_geometry_params(theta_i, theta_s, phi_i, phi_s, k);
    
    % Build quadrature grid
    quad = build_quadrature(kl, n_points, nmax);
    
    % Pre-compute constants
    constants = precompute_constants(correlation_type, nmax);
    
    % Create roughness spectrum provider
    wn_provider = make_Wn_provider(sigma, kl, correlation_type, constants);
    
    % Compute vertical wavenumbers
    U = quad.U;
    V = quad.V;
    q1 = sqrt(max(k^2 - (U.^2 + V.^2), 0));
    q2 = sqrt(eps_r * k^2 - (U.^2 + V.^2));
    
    % Radiation condition mask
    qmin = 1e-6;
    rad = (real(q1) > qmin) | (real(q2) > qmin);
    
    % Quadrature weights
    W2D = quad.wu' * quad.wv;
    
    % Build propagators
    propagators = build_propagators(U, V, q1, q2, k, eps_r, geom, polarization);
    
    % Build Kirchhoff-complementary terms
    K1 = build_gkc1(U, V, geom, q1, sigma, kl, wn_provider, nmax, constants);
    K2 = build_gkc2(U, V, geom, q1, sigma, kl, wn_provider, nmax, constants);
    K3 = build_gkc3(U, V, geom, q1, sigma, kl, wn_provider, nmax, constants);
    
    % Build complementary terms
    C1 = build_gc_block1(U, V, geom, q1, q1, sigma, kl, wn_provider, nmax, constants);
    C2 = build_gc_block2(U, V, geom, q1, q1, sigma, kl, wn_provider, nmax, constants);
    
    % Assemble integrands
    P = propagators;
    
    % Kirchhoff-complementary integrand
    Int_kc = abs(P.Fp).^2 .* K1 + abs(P.Fm).^2 .* K2 + abs(P.Gp).^2 .* K3;
    
    % Complementary integrand (all 14 terms)
    Int_c = zeros(size(U));
    Int_c = Int_c + abs(P.Fp).^2 .* C1.gc1;
    Int_c = Int_c + (P.Fp .* conj(P.Fm)) .* C1.gc2;
    Int_c = Int_c + (P.Fm .* conj(P.Fp)) .* C1.gc3;
    Int_c = Int_c + abs(P.Fm).^2 .* C1.gc4;
    Int_c = Int_c + abs(P.Gp).^2 .* C1.gc5;
    Int_c = Int_c + (P.Gp .* conj(P.Gm)) .* C1.gc6;
    Int_c = Int_c + (P.Gm .* conj(P.Gp)) .* C1.gc7;
    Int_c = Int_c + abs(P.Gm).^2 .* C1.gc8;
    Int_c = Int_c + (P.Fp .* conj(P.Gp)) .* C2.gc9;
    Int_c = Int_c + (P.Fp .* conj(P.Gm)) .* C2.gc10;
    Int_c = Int_c + (P.Fm .* conj(P.Gp)) .* C2.gc11;
    Int_c = Int_c + (P.Fm .* conj(P.Gm)) .* C2.gc12;
    Int_c = Int_c + (P.Gp .* conj(P.Fp)) .* C2.gc13;
    Int_c = Int_c + (P.Gm .* conj(P.Fp)) .* C2.gc14;
    
    % Apply radiation condition and take real part
    % NOTE: Do NOT take abs() - allow negative values to cancel during integration
    Ikc_real = real(Int_kc) .* rad;
    Ic_real = real(Int_c) .* rad;
    
    % Integrate
    val_kc = sum(sum(Ikc_real .* W2D));
    val_c = sum(sum(Ic_real .* W2D));
    
    % Multiple scattering coefficient
    sigma_ms = (k^2 / (8.0 * pi)) * val_kc + (k^2 / (64.0 * pi)) * val_c;
    sigma_ms = max(real(sigma_ms), 0.0);
end


% ============================================================================
% GEOMETRY PREPARATION
% ============================================================================

function geom = prepare_geometry_params(theta_i, theta_s, phi_i, phi_s, k)
    % Prepare geometry parameters with precomputed trigonometric values
    
    geom.theta_i = theta_i;
    geom.theta_s = theta_s;
    geom.phi_i = phi_i;
    geom.phi_s = phi_s;
    
    geom.sin_theta_i = sin(theta_i);
    geom.cos_theta_i = cos(theta_i);
    geom.sin_theta_s = sin(theta_s);
    geom.cos_theta_s = cos(theta_s);
    geom.sin_phi_i = sin(phi_i);
    geom.cos_phi_i = cos(phi_i);
    geom.sin_phi_s = sin(phi_s);
    geom.cos_phi_s = cos(phi_s);
    
    geom.kx = k * geom.sin_theta_i * geom.cos_phi_i;
    geom.ky = k * geom.sin_theta_i * geom.sin_phi_i;
    geom.kz = k * geom.cos_theta_i;
    geom.ksx = k * geom.sin_theta_s * geom.cos_phi_s;
    geom.ksy = k * geom.sin_theta_s * geom.sin_phi_s;
    geom.ksz = k * geom.cos_theta_s;
end


% ============================================================================
% QUADRATURE GRID
% ============================================================================

function quad = build_quadrature(kl, n_points, nmax)
    % Build quadrature grid for spectral integration
    
    umax = 10.0 / max(kl, 1e-6);
    grid = linspace(-umax, umax, n_points);
    [U, V] = meshgrid(grid, grid);
    
    if length(grid) > 1
        du = grid(2) - grid(1);
    else
        du = 0.0;
    end
    dv = du;
    
    wu = ones(size(grid)) * du;
    wv = ones(size(grid)) * dv;
    
    if length(grid) > 1
        wu(1) = 0.5 * du;
        wu(end) = 0.5 * du;
        wv(1) = 0.5 * dv;
        wv(end) = 0.5 * dv;
    end
    
    quad.U = U;
    quad.V = V;
    quad.wu = wu;
    quad.wv = wv;
    quad.Nmax = nmax;
end


% ============================================================================
% CONSTANTS PRE-COMPUTATION
% ============================================================================

function constants = precompute_constants(correlation_type, nmax)
    % Pre-compute constants for acceleration
    %
    % CRITICAL NORMALIZATION FACTOR: (2π)^10
    % This factor is essential for correct magnitude scaling of the
    % exponential correlation spectrum in multiple scattering.
    
    % Pre-compute factorials for series summation
    constants.factorials = zeros(nmax, 1);
    fact = 1.0;
    for n = 1:nmax
        fact = fact * n;
        constants.factorials(n) = fact;
    end
    
    % Pre-compute normalization factor for exponential correlation
    % This is the critical (2π)^10 factor that was missing
    if correlation_type == 2  % Exponential
        constants.two_pi_power = (2.0 * pi)^10;
    else
        constants.two_pi_power = 1.0;  % Not used for Gaussian
    end
end

% ============================================================================
% ROUGHNESS SPECTRUM PROVIDER
% ============================================================================

function wn_provider = make_Wn_provider(sigma, kl, correlation_type, constants)
    % Create roughness spectrum provider
    
    sigma2 = sigma^2;
    
    if correlation_type == 1  % Gaussian
        wn_provider = @(u, v, n) gaussian_spectrum(u, v, n, sigma2, kl);
    else  % Exponential
        two_pi_power = constants.two_pi_power;
        wn_provider = @(u, v, n) exponential_spectrum(u, v, n, sigma2, kl, two_pi_power);
    end
end

function W = gaussian_spectrum(u, v, n, sigma2, kl)
    % Gaussian correlation spectrum
    factor = kl^2 / max(n, 1);
    exp_arg = -(kl^2 / (4.0 * max(n, 1))) * (u.^2 + v.^2);
    W = sigma2 * (factor / (4.0 * pi)) * exp(exp_arg);
end

function W = exponential_spectrum(u, v, n, sigma2, kl, two_pi_power)
    % Exponential correlation spectrum
    denom = 1.0 + ((kl * sqrt(u.^2 + v.^2)) / max(n, 1)).^2;
    W = two_pi_power * sigma2 * (kl / max(n, 1))^2 * denom.^(-1.5);
end


% ============================================================================
% PROPAGATOR COMPUTATION
% ============================================================================

function propagators = build_propagators(U, V, q1, q2, k, eps_r, geom, pol)
    % Build upward and downward propagators for given polarization
    
    % Spectral angles
    [cos_phi, sin_phi] = spectral_angles(U, V);
    
    % Compute C and B coefficients
    C_air = compute_C_coeffs(q1, geom, cos_phi, sin_phi, U, V);
    C_soil = compute_C_coeffs(q2, geom, cos_phi, sin_phi, U, V);
    B_air = compute_B_coeffs(q1, geom, cos_phi, sin_phi, U, V);
    B_soil = compute_B_coeffs(q2, geom, cos_phi, sin_phi, U, V);
    
    % Fresnel coefficients
    [Rh, Rv] = fresnel_coeffs(eps_r, q1, q2);
    R = 0.5 * (Rv - Rh);
    
    % Material parameters
    mu_r = 1.0;
    u_r = 1.0;
    inv_q1 = safe_inverse(q1);
    inv_q2 = safe_inverse(q2);
    
    % Compute upward propagators
    pol = lower(pol);
    if strcmp(pol, 'vv')
        coeff = C_air;
        coeff_t = C_soil;
        Fp_plus = ...
            -(1 - Rv) .* (1 + Rv) .* inv_q1 .* coeff.C1 + ...
            (1 - Rv) .* (1 - Rv) .* inv_q1 .* coeff.C2 + ...
            (1 - Rv) .* (1 + Rv) .* inv_q1 .* coeff.C3 + ...
            (1 + Rv) .* (1 - Rv) .* inv_q1 .* coeff.C4 + ...
            (1 + Rv) .* (1 + Rv) .* inv_q1 .* coeff.C5 + ...
            (1 + Rv) .* (1 - Rv) .* inv_q1 .* coeff.C6;
        Gp_plus = ...
            (1 + Rv) .* (1 + Rv) * u_r .* inv_q2 .* coeff_t.C1 - ...
            (1 + Rv) .* (1 - Rv) .* inv_q2 .* coeff_t.C2 - ...
            (1 + Rv * eps_r) .* (1 + Rv) .* inv_q2 .* coeff_t.C3 - ...
            (1 - Rv) * eps_r .* (1 - Rv) .* inv_q2 .* coeff_t.C4 - ...
            (1 - Rv) .* (1 + Rv) .* inv_q2 .* coeff_t.C5 - ...
            (1 - Rv) .* (1 - Rv) * u_r .* inv_q2 .* coeff_t.C6;
    elseif strcmp(pol, 'hh')
        coeff = C_air;
        coeff_t = C_soil;
        Fp_plus = ...
            (1 - Rh) .* (1 + Rh) .* inv_q1 .* coeff.C1 - ...
            (1 - Rh) .* (1 - Rh) .* inv_q1 .* coeff.C2 - ...
            (1 - Rh) .* (1 + Rh) .* inv_q1 .* coeff.C3 - ...
            (1 + Rh) .* (1 - Rh) .* inv_q1 .* coeff.C4 - ...
            (1 + Rh) .* (1 + Rh) .* inv_q1 .* coeff.C5 - ...
            (1 + Rh) .* (1 - Rh) .* inv_q1 .* coeff.C6;
        Gp_plus = ...
            -(1 + Rh) .* (1 + Rh) * eps_r .* inv_q2 .* coeff_t.C1 + ...
            (1 + Rh) .* (1 - Rh) .* inv_q2 .* coeff_t.C2 + ...
            (1 + Rh) .* (1 + Rh) * u_r .* inv_q2 .* coeff_t.C3 + ...
            (1 - Rh) .* (1 - Rh) * u_r .* inv_q2 .* coeff_t.C4 + ...
            (1 - Rh) .* (1 + Rh) .* inv_q2 .* coeff_t.C5 + ...
            (1 - Rh) .* (1 - Rh) * eps_r .* inv_q2 .* coeff_t.C6;
    else  % cross-pol
        coeff = B_air;
        coeff_t = B_soil;
        Fp_plus = ...
            (1 - R) .* (1 + R) .* inv_q1 .* coeff.B1 - ...
            (1 - R) .* (1 - R) .* inv_q1 .* coeff.B2 - ...
            (1 - R) .* (1 + R) .* inv_q1 .* coeff.B3 + ...
            (1 + R) .* (1 - R) .* inv_q1 .* coeff.B4 + ...
            (1 + R) .* (1 + R) .* inv_q1 .* coeff.B5 + ...
            (1 + R) .* (1 - R) .* inv_q1 .* coeff.B6;
        Gp_plus = ...
            -(1 + R) .* (1 + R) * mu_r .* inv_q2 .* coeff_t.B1 + ...
            (1 + R) .* (1 - R) .* inv_q2 .* coeff_t.B2 + ...
            (1 + R * eps_r) .* (1 + R) .* inv_q2 .* coeff_t.B3 - ...
            (1 - R) * eps_r .* (1 - R) .* inv_q2 .* coeff_t.B4 - ...
            (1 - R) .* (1 + R) .* inv_q2 .* coeff_t.B5 - ...
            (1 - R) .* (1 - R) * mu_r .* inv_q2 .* coeff_t.B6;
    end
    
    % Compute downward propagators
    [Fm_plus, Gm_plus] = compute_downward_propagators(U, V, q1, q2, k, eps_r, geom, pol, cos_phi, sin_phi);
    
    propagators.Fp = Fp_plus;
    propagators.Fm = Fm_plus;
    propagators.Gp = Gp_plus;
    propagators.Gm = Gm_plus;
end

function [Fm, Gm] = compute_downward_propagators(U, V, q1, q2, k, eps_r, geom, pol, cos_phi, sin_phi)
    % Compute downward propagators
    
    % Use negative q for downward propagation
    [Rh, Rv] = fresnel_coeffs(eps_r, -q1, -q2);
    R = 0.5 * (Rv - Rh);
    mu_r = 1.0;
    u_r = 1.0;
    inv_q1 = safe_inverse(-q1);
    inv_q2 = safe_inverse(-q2);
    
    C_air = compute_C_coeffs(-q1, geom, cos_phi, sin_phi, U, V);
    C_soil = compute_C_coeffs(-q2, geom, cos_phi, sin_phi, U, V);
    B_air = compute_B_coeffs(-q1, geom, cos_phi, sin_phi, U, V);
    B_soil = compute_B_coeffs(-q2, geom, cos_phi, sin_phi, U, V);
    
    if strcmp(pol, 'vv')
        coeff = C_air;
        coeff_t = C_soil;
        Fm = ...
            -(1 - Rv) .* (1 + Rv) .* inv_q1 .* coeff.C1 + ...
            (1 - Rv) .* (1 - Rv) .* inv_q1 .* coeff.C2 + ...
            (1 - Rv) .* (1 + Rv) .* inv_q1 .* coeff.C3 + ...
            (1 + Rv) .* (1 - Rv) .* inv_q1 .* coeff.C4 + ...
            (1 + Rv) .* (1 + Rv) .* inv_q1 .* coeff.C5 + ...
            (1 + Rv) .* (1 - Rv) .* inv_q1 .* coeff.C6;
        Gm = ...
            (1 + Rv) .* (1 + Rv) * u_r .* inv_q2 .* coeff_t.C1 - ...
            (1 + Rv) .* (1 - Rv) .* inv_q2 .* coeff_t.C2 - ...
            (1 + Rv * eps_r) .* (1 + Rv) .* inv_q2 .* coeff_t.C3 - ...
            (1 - Rv) * eps_r .* (1 - Rv) .* inv_q2 .* coeff_t.C4 - ...
            (1 - Rv) .* (1 + Rv) .* inv_q2 .* coeff_t.C5 - ...
            (1 - Rv) .* (1 - Rv) * u_r .* inv_q2 .* coeff_t.C6;
    elseif strcmp(pol, 'hh')
        coeff = C_air;
        coeff_t = C_soil;
        Fm = ...
            (1 - Rh) .* (1 + Rh) .* inv_q1 .* coeff.C1 - ...
            (1 - Rh) .* (1 - Rh) .* inv_q1 .* coeff.C2 - ...
            (1 - Rh) .* (1 + Rh) .* inv_q1 .* coeff.C3 - ...
            (1 + Rh) .* (1 - Rh) .* inv_q1 .* coeff.C4 - ...
            (1 + Rh) .* (1 + Rh) .* inv_q1 .* coeff.C5 - ...
            (1 + Rh) .* (1 - Rh) .* inv_q1 .* coeff.C6;
        Gm = ...
            -(1 + Rh) .* (1 + Rh) * eps_r .* inv_q2 .* coeff_t.C1 + ...
            (1 + Rh) .* (1 - Rh) .* inv_q2 .* coeff_t.C2 + ...
            (1 + Rh) .* (1 + Rh) .* inv_q2 .* coeff_t.C3 + ...
            (1 - Rh) .* (1 - Rh) .* inv_q2 .* coeff_t.C4 + ...
            (1 - Rh) .* (1 + Rh) .* inv_q2 .* coeff_t.C5 + ...
            (1 - Rh) .* (1 - Rh) * eps_r .* inv_q2 .* coeff_t.C6;
    else  % cross-pol
        coeff = B_air;
        coeff_t = B_soil;
        Fm = ...
            (1 - R) .* (1 + R) .* inv_q1 .* coeff.B1 - ...
            (1 - R) .* (1 - R) .* inv_q1 .* coeff.B2 - ...
            (1 - R) .* (1 + R) .* inv_q1 .* coeff.B3 + ...
            (1 + R) .* (1 - R) .* inv_q1 .* coeff.B4 + ...
            (1 + R) .* (1 + R) .* inv_q1 .* coeff.B5 + ...
            (1 + R) .* (1 - R) .* inv_q1 .* coeff.B6;
        Gm = ...
            -(1 + R) .* (1 + R) * mu_r .* inv_q2 .* coeff_t.B1 + ...
            (1 + R) .* (1 - R) .* inv_q2 .* coeff_t.B2 + ...
            (1 + R * eps_r) .* (1 + R) .* inv_q2 .* coeff_t.B3 - ...
            (1 - R) * eps_r .* (1 - R) .* inv_q2 .* coeff_t.B4 - ...
            (1 - R) .* (1 + R) .* inv_q2 .* coeff_t.B5 - ...
            (1 - R) .* (1 - R) * mu_r .* inv_q2 .* coeff_t.B6;
    end
end


% ============================================================================
% C AND B COEFFICIENTS (Appendix C of Yang et al. 2017)
% ============================================================================

function [cos_phi, sin_phi] = spectral_angles(U, V)
    % Compute spectral domain angles
    rho = sqrt(U.^2 + V.^2);
    cos_phi = ones(size(U));
    sin_phi = zeros(size(V));
    mask = rho > 0;
    cos_phi(mask) = U(mask) ./ rho(mask);
    sin_phi(mask) = V(mask) ./ rho(mask);
end

function C = compute_C_coeffs(q, geom, cos_phi, sin_phi, U, V)
    % Compute C coefficients for VV/HH polarizations (Eqs C1-C6)
    
    cos_phi_s = geom.cos_phi_s;
    sin_phi_s = geom.sin_phi_s;
    cos_theta = geom.cos_theta_s;
    sin_theta = geom.sin_theta_s;
    z_x = geom.sin_theta_i * geom.cos_phi_i;
    z_y = geom.sin_theta_i * geom.sin_phi_i;
    zp_x = geom.sin_theta_s * geom.cos_phi_s;
    zp_y = geom.sin_theta_s * geom.sin_phi_s;
    
    C.C1 = -cos_phi_s * (-cos_phi - z_x * zp_x * cos_phi - z_x * zp_y * sin_phi) + ...
           sin_phi_s * (sin_phi + zp_x * z_y * cos_phi + z_y * zp_y * sin_phi);
    
    C.C2 = -cos_phi_s * (-q * cos_theta * cos_phi - U * z_x * cos_theta - ...
           V * zp_y * cos_theta * cos_phi - q * zp_x * sin_theta - ...
           U * z_x * zp_x * sin_theta - V * z_x * zp_y * sin_theta - ...
           V * z_x * cos_theta * sin_phi + V * zp_x * cos_theta * sin_phi) + ...
           sin_phi_s * (U * z_y * cos_theta * cos_phi - U * zp_y * cos_theta + ...
           U * zp_x * z_y * sin_theta + q * zp_y * sin_theta + ...
           V * z_y * zp_y * sin_theta + q * cos_theta * sin_phi + ...
           U * zp_x * cos_theta * sin_phi + V * z_y * cos_theta * sin_phi);
    
    C.C3 = cos_phi_s * (U * zp_x * cos_theta * cos_phi - q * z_x * zp_x * cos_theta * cos_phi - ...
           U * sin_theta + q * z_x * sin_theta + U * zp_y * cos_theta * sin_phi - ...
           q * z_x * zp_y * cos_theta * sin_phi) + ...
           sin_phi_s * (V * zp_x * cos_theta * cos_phi - q * zp_x * z_y * cos_theta * cos_phi - ...
           V * sin_theta + q * z_y * sin_theta + V * zp_y * cos_theta * sin_phi - ...
           q * z_y * zp_y * cos_theta * sin_phi);
    
    C.C4 = sin_theta * (-z_x * cos_theta * cos_phi - z_x * zp_x * sin_theta - ...
           z_y * zp_y * sin_theta - z_y * cos_theta * sin_phi) - ...
           cos_theta * cos_phi_s * (-cos_theta * cos_phi - z_y * zp_y * cos_theta * cos_phi - ...
           zp_x * sin_theta + zp_x * z_y * cos_theta * sin_phi) - ...
           cos_theta * sin_phi_s * (z_x * zp_y * cos_theta * cos_phi - zp_y * sin_theta - ...
           cos_theta * sin_phi - z_x * zp_x * cos_theta * sin_phi);
    
    C.C5 = sin_theta * (q * z_x * cos_phi + U * z_x * zp_x * cos_phi + ...
           V * zp_x * z_y * cos_phi + q * z_y * sin_phi + U * z_x * zp_y * sin_phi + ...
           V * z_y * zp_y * sin_phi) - ...
           cos_theta * cos_phi_s * (q * cos_phi + U * zp_x * cos_phi + V * z_y * cos_phi - ...
           U * z_y * sin_phi + U * zp_y * sin_phi) - ...
           cos_theta * sin_phi_s * (-V * z_x * cos_phi + V * zp_x * cos_phi + ...
           q * sin_phi + U * z_x * sin_phi + V * zp_y * sin_phi);
    
    C.C6 = sin_theta * (V * z_x * zp_y * cos_phi - U * z_y * zp_y * cos_phi - ...
           V * z_x * zp_x * sin_phi + U * zp_x * z_y * sin_phi) + ...
           cos_theta * cos_phi_s * (-V * zp_y * cos_phi + q * z_y * zp_y * cos_phi + ...
           V * zp_x * sin_phi - q * zp_x * z_y * sin_phi) + ...
           cos_theta * sin_phi_s * (U * zp_y * cos_phi - q * z_x * zp_y * cos_phi - ...
           U * zp_x * sin_phi + q * z_x * zp_x * sin_phi);
end

function B = compute_B_coeffs(q, geom, cos_phi, sin_phi, U, V)
    % Compute B coefficients for HV/VH polarizations (Eqs C7-C12)
    
    cos_phi_s = geom.cos_phi_s;
    sin_phi_s = geom.sin_phi_s;
    cos_theta = geom.cos_theta_s;
    sin_theta = geom.sin_theta_s;
    z_x = geom.sin_theta_i * geom.cos_phi_i;
    z_y = geom.sin_theta_i * geom.sin_phi_i;
    zp_x = geom.sin_theta_s * geom.cos_phi_s;
    zp_y = geom.sin_theta_s * geom.sin_phi_s;
    
    B.B1 = sin_theta * (-z_y * cos_phi + z_x * sin_phi) - ...
           cos_theta * cos_phi_s * (zp_x * z_y * cos_phi + sin_phi + z_y * zp_y * sin_phi) - ...
           cos_theta * sin_phi_s * (-cos_phi - z_x * zp_x * cos_phi - z_x * zp_y * sin_phi);
    
    B.B2 = sin_theta * (-q * z_y * cos_theta * cos_phi - U * z_x * zp_y * cos_theta * cos_phi - ...
           V * z_y * zp_y * cos_theta * cos_phi - q * zp_x * z_y * sin_theta + ...
           q * z_x * zp_y * sin_theta + q * z_x * cos_theta * sin_phi + ...
           U * z_x * zp_x * cos_theta * sin_phi + V * zp_x * z_y * cos_theta * sin_phi) - ...
           cos_theta * cos_phi_s * (U * z_y * cos_theta * cos_phi - U * zp_y * cos_theta * cos_phi + ...
           U * zp_x * z_y * sin_theta + q * zp_y * sin_theta + V * z_y * zp_y * sin_theta + ...
           q * cos_theta * sin_phi + U * zp_x * cos_theta * sin_phi + V * z_y * cos_theta * sin_phi) - ...
           cos_theta * sin_phi_s * (-q * cos_theta * cos_phi - U * z_x * cos_theta * cos_phi - ...
           V * zp_y * cos_theta * cos_phi - q * zp_x * sin_theta - U * z_x * zp_x * sin_theta - ...
           V * z_x * zp_y * sin_theta - V * z_x * cos_theta * sin_phi + V * zp_x * cos_theta * sin_phi);
    
    B.B3 = sin_theta * (V * z_x * zp_x * cos_theta * cos_phi - U * z_y * zp_x * cos_theta * cos_phi - ...
           V * z_x * sin_theta + U * z_y * sin_theta + V * z_x * zp_y * cos_theta * sin_phi - ...
           U * z_y * zp_y * cos_theta * sin_phi) + ...
           cos_theta * cos_phi_s * (-V * zp_x * cos_theta * cos_phi + q * zp_x * z_y * cos_theta * cos_phi + ...
           V * sin_theta - q * z_y * sin_theta - V * zp_y * cos_theta * sin_phi + ...
           q * z_y * zp_y * cos_theta * sin_phi) + ...
           cos_theta * sin_phi_s * (U * zp_x * cos_theta * cos_phi - q * z_x * zp_x * cos_theta * cos_phi - ...
           U * sin_theta + q * z_x * sin_theta + U * zp_y * cos_theta * sin_phi - ...
           q * z_x * zp_y * cos_theta * sin_phi);
    
    B.B4 = -cos_phi_s * (z_x * zp_y * cos_theta * cos_phi - zp_y * sin_theta - ...
           cos_theta * sin_phi - z_x * zp_x * cos_theta * sin_phi) + ...
           sin_phi_s * (-cos_theta * cos_phi - z_y * zp_y * cos_theta * cos_phi - ...
           zp_x * sin_theta + zp_x * z_y * cos_theta * sin_phi);
    
    B.B5 = -cos_phi_s * (-V * z_x * cos_phi + V * zp_x * cos_phi + q * sin_phi + ...
           U * z_x * sin_phi - V * zp_y * sin_phi) + ...
           sin_phi_s * (q * cos_phi + U * zp_x * cos_phi + V * z_y * cos_phi - ...
           U * z_y * sin_phi + U * zp_y * sin_phi);
    
    B.B6 = cos_phi_s * (U * zp_y * cos_phi - q * z_x * zp_y * cos_phi - ...
           U * zp_x * sin_phi + q * z_x * zp_x * sin_phi) + ...
           sin_phi_s * (V * zp_y * cos_phi - q * z_y * zp_y * cos_phi - ...
           V * zp_x * sin_phi + q * zp_x * z_y * sin_phi);
end


% ============================================================================
% UTILITY FUNCTIONS
% ============================================================================

function inv_x = safe_inverse(x)
    % Safe inverse avoiding singularities
    safe_eps = 1e-8;
    mask = abs(x) < safe_eps;
    x_safe = x + mask * (safe_eps + 1i * safe_eps);
    inv_x = 1.0 ./ x_safe;
end

function [Rh, Rv] = fresnel_coeffs(eps_r, q1, q2)
    % Compute Fresnel reflection coefficients
    Rh = safe_div(q1 - q2, q1 + q2);
    Rv = safe_div(eps_r * q1 - q2, eps_r * q1 + q2);
end

function result = safe_div(a, b)
    % Safe division
    safe_eps = 1e-8;
    mask = abs(b) < safe_eps;
    b_safe = b + mask * (safe_eps + 1i * safe_eps);
    result = a ./ b_safe;
end


% ============================================================================
% KIRCHHOFF-COMPLEMENTARY TERMS (Appendix A, Eqs A1-A3)
% ============================================================================

function K1 = build_gkc1(U, V, geom, q, sigma, kl, wn_provider, Nmax, constants)
    % Build Kirchhoff-complementary term K1 (Eq A1)
    
    sigma2 = sigma^2;
    kz = geom.kz;
    ksz = geom.ksz;
    kx = geom.kx;
    ky = geom.ky;
    ksx = geom.ksx;
    ksy = geom.ksy;
    
    expo = exp(-sigma2 * (ksz^2 + kz^2 + ksz * kz + q.^2 - q * ksz + q * kz));
    a_m = sigma2 * (kz + q) * (ksz + kz);
    a_n = sigma2 * (ksz - q) * (ksz + kz);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    K1 = expo .* sum_m .* sum_n;
end

function K2 = build_gkc2(U, V, geom, q, sigma, kl, wn_provider, Nmax, constants)
    % Build Kirchhoff-complementary term K2 (Eq A2)
    
    sigma2 = sigma^2;
    kz = geom.kz;
    ksz = geom.ksz;
    kx = geom.kx;
    ky = geom.ky;
    ksx = geom.ksx;
    ksy = geom.ksy;
    
    expo = exp(-sigma2 * (ksz^2 + kz^2 + ksz * kz + q.^2 - q * ksz + q * kz));
    a_m = sigma2 * (kz + q) * (ksz + kz);
    a_n = -sigma2 * (ksz - q) * (kz + q);
    sum_m = series_sum(a_m, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    K2 = expo .* sum_m .* sum_n;
end

function K3 = build_gkc3(U, V, geom, q, sigma, kl, wn_provider, Nmax, constants)
    % Build Kirchhoff-complementary term K3 (Eq A3)
    
    sigma2 = sigma^2;
    kz = geom.kz;
    ksz = geom.ksz;
    kx = geom.kx;
    ky = geom.ky;
    ksx = geom.ksx;
    ksy = geom.ksy;
    
    expo = exp(-sigma2 * (ksz^2 + kz^2 + ksz * kz + q.^2 - q * ksz + q * kz));
    a_m = -sigma2 * (ksz - q) * (kz + q);
    a_n = sigma2 * (ksz - q) * (ksz + kz);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    K3 = expo .* sum_m .* sum_n;
end


% ============================================================================
% COMPLEMENTARY TERMS BLOCK 1 (Appendix A, Eqs A4-A11)
% ============================================================================

function C1 = build_gc_block1(U, V, geom, q, qp, sigma, kl, wn_provider, Nmax, constants)
    % Build complementary block 1 (gc1-gc8)
    
    sigma2 = sigma^2;
    kz = geom.kz;
    ksz = geom.ksz;
    kx = geom.kx;
    ky = geom.ky;
    ksx = geom.ksx;
    ksy = geom.ksy;
    
    expo_func = @(q_, qp_) exp(-sigma2 * (ksz^2 + kz^2 + q_.^2 + qp_.^2 - (ksz - kz) * (q_ + qp_)));
    
    % gc1 (A4)
    a_n = sigma2 * (ksz - q) * (ksz - qp);
    a_m = sigma2 * (kz + q) * (kz + qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc1 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc2 (A5)
    a_n = sigma2 * (ksz - q) * (kz + qp);
    a_m = sigma2 * (kz + q) * (ksz - qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc2 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc3 (A6)
    a_n = sigma2 * (ksz - q) * (kz + qp);
    a_m = sigma2 * (kz + q) * (kz + qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc3 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc4 (A7)
    a_n = sigma2 * (ksz - q) * (kz + qp);
    a_m = -sigma2 * (ksz - q) * (kz + q);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc4 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc5 (A8)
    a_n = sigma2 * (kz + q) * (kz + qp);
    a_m = -sigma2 * (ksz - q) * (kz + q);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax, constants);
    C1.gc5 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc6 (A9)
    a_n = sigma2 * (ksz - q) * (ksz - qp);
    a_m = sigma2 * (kz + q) * (ksz - qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc6 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc7 (A10)
    a_n = sigma2 * (ksz - q) * (ksz - qp);
    a_m = -sigma2 * (ksz - q) * (kz + q);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C1.gc7 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc8 (A11)
    a_n = sigma2 * (kz + q) * (ksz - qp);
    a_m = -sigma2 * (ksz - q) * (kz + q);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax, constants);
    C1.gc8 = expo_func(q, qp) .* sum_n .* sum_m;
end


% ============================================================================
% COMPLEMENTARY TERMS BLOCK 2 (Appendix A, Eqs A12-A17)
% ============================================================================

function C2 = build_gc_block2(U, V, geom, q, qp, sigma, kl, wn_provider, Nmax, constants)
    % Build complementary block 2 (gc9-gc14)
    
    sigma2 = sigma^2;
    kz = geom.kz;
    ksz = geom.ksz;
    kx = geom.kx;
    ky = geom.ky;
    ksx = geom.ksx;
    ksy = geom.ksy;
    
    expo_func = @(q_, qp_) exp(-sigma2 * (ksz^2 + kz^2 + q_.^2 + qp_.^2 - (ksz - kz) * (q_ + qp_)));
    
    % gc9 (A12)
    a_n = sigma2 * (kz + q) * (ksz - qp);
    a_m = sigma2 * (kz + q) * (kz + qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C2.gc9 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc10 (A13)
    a_n = sigma2 * (kz + q) * (ksz - qp);
    a_m = -sigma2 * (ksz - qp) * (kz + qp);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C2.gc10 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc11 (A14)
    a_n = sigma2 * (kz + q) * (kz + qp);
    a_m = -sigma2 * (ksz - qp) * (kz + qp);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax, constants);
    C2.gc11 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc12 (A15)
    a_n = sigma2 * (ksz - q) * (ksz - qp);
    a_m = sigma2 * (ksz - q) * (kz + qp);
    sum_n = series_sum(a_n, ksx + U, ksy + V, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C2.gc12 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc13 (A16)
    a_n = sigma2 * (ksz - q) * (ksz - qp);
    a_m = -sigma2 * (ksz - qp) * (kz + qp);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, kx + U, ky + V, wn_provider, Nmax, constants);
    C2.gc13 = expo_func(q, qp) .* sum_n .* sum_m;
    
    % gc14 (A17)
    a_n = sigma2 * (ksz - q) * (kz + qp);
    a_m = -sigma2 * (ksz - qp) * (kz + qp);
    sum_n = series_sum(a_n, kx - ksx, ky - ksy, wn_provider, Nmax, constants);
    sum_m = series_sum(a_m, ksx + U, ksy + V, wn_provider, Nmax, constants);
    C2.gc14 = expo_func(q, qp) .* sum_n .* sum_m;
end


% ============================================================================
% SERIES SUMMATION
% ============================================================================

function result = series_sum(coeff, arg_x, arg_y, wn_provider, Nmax, constants)
    % Compute series summation with roughness spectrum
    
    result = zeros(size(arg_x));
    factorials = constants.factorials;
    
    for n = 1:Nmax
        Wn = wn_provider(arg_x, arg_y, n);
        factorial_n = factorials(n);
        result = result + (coeff.^n / factorial_n) .* Wn;
    end
end
