function [Rvi, Rhi, Rvhi, Rvl, Rhl, Rvhl, rv0, rh0] = fresnel_coefficients(eps_r, theta_i, theta_s, phi_s)
% FRESNEL_COEFFICIENTS Compute Fresnel reflection coefficients with bug fixes
%
% Computes Fresnel coefficients at incident angle, specular angle, and normal incidence
% with proper branch selection for lossy media.
%
% Inputs:
%   eps_r   - Complex relative permittivity
%   theta_i - Incident angle (radians)
%   theta_s - Scattered angle (radians)
%   phi_s   - Scattered azimuth angle (radians)
%
% Outputs:
%   Rvi, Rhi, Rvhi - Reflection coefficients at incident angle
%   Rvl, Rhl, Rvhl - Reflection coefficients at specular angle
%   rv0, rh0       - Reflection coefficients at normal incidence
%
% Bug fixes applied:
%   1. Proper branch selection for lossy media (Im(stem) >= 0)
%   2. Correct normal incidence: rh0 = rv0 (not -rv0)

    mu_r = 1.0;  % Non-magnetic soil
    
    % Trigonometric quantities
    si = sin(theta_i);
    cs = cos(theta_i);
    si2 = si * si;
    sis = sin(theta_s);
    css = cos(theta_s);
    csfs = cos(phi_s);
    
    %% Fresnel coefficients at incident angle
    % BUG FIX 1: Proper branch selection for lossy media
    stem = sqrt(eps_r * mu_r - si2);
    if imag(stem) < 0
        stem = -stem;  % Ensure Im(stem) >= 0 for decaying wave
    end
    
    % Vertical (TM) and Horizontal (TE) polarization
    Rvi = (eps_r * cs - stem) / (eps_r * cs + stem);
    Rhi = (mu_r * cs - stem) / (mu_r * cs + stem);
    Rvhi = (Rvi - Rhi) / 2.0;
    
    %% Fresnel coefficients at specular angle
    % Compute specular half-angle
    csl = sqrt((1.0 + cs * css - si * sis * csfs) / 2.0);
    sil = sqrt(1.0 - csl * csl);
    
    % BUG FIX 1: Proper branch selection for lossy media
    steml = sqrt(eps_r * mu_r - sil * sil);
    if imag(steml) < 0
        steml = -steml;  % Ensure Im(steml) >= 0 for decaying wave
    end
    
    Rvl = (eps_r * csl - steml) / (eps_r * csl + steml);
    Rhl = (mu_r * csl - steml) / (mu_r * csl + steml);
    Rvhl = (Rvl - Rhl) / 2.0;
    
    %% Normal incidence coefficients
    % BUG FIX 2: Both polarizations use same value at normal incidence
    sqrt_er = sqrt(eps_r);
    r0 = (sqrt_er - 1.0) / (sqrt_er + 1.0);
    rv0 = r0;
    rh0 = r0;  % CORRECTED: was -r0 in original MATLAB
end
