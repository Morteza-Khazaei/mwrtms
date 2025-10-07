function [Tfv, Tfh] = transition_function(eps_r, theta_i, ks, cs, spectra, n_terms)
% TRANSITION_FUNCTION Compute AIEM transition function with bug fixes
%
% Computes transition function for V and H polarizations using the
% Wu & Fung (1992) formulation with corrections.
%
% Inputs:
%   eps_r   - Complex relative permittivity
%   theta_i - Incident angle (radians)
%   ks      - Normalized RMS height (k * sigma)
%   cs      - cos(theta_i)
%   spectra - Roughness spectrum values W^(n) [1 x n_terms]
%   n_terms - Number of spectral terms
%
% Outputs:
%   Tfv - Transition function for V polarization
%   Tfh - Transition function for H polarization
%
% Bug fixes applied:
%   1. rh0 = rv0 (not -rv0) at normal incidence
%   2. Use rh0 (not rv0) in H-polarization path
%
% Note: This is the LEGACY Wu & Fung method. For <1 dB RMSE vs NMM3D,
%       implement the new S_p/S_p^(0) method from the bug report.

    si = sin(theta_i);
    s2 = si * si;
    
    % BUG FIX 1: Normal incidence coefficients (both same)
    sqrt_er = sqrt(eps_r);
    rv0 = (sqrt_er - 1.0) / (sqrt_er + 1.0);
    rh0 = rv0;  % CORRECTED: was -rv0 in original MATLAB
    
    % Transmitted wave vector component
    rt = sqrt(eps_r - s2);
    
    % Transition factors
    denom_v = cs * rt;
    if abs(denom_v) < 1e-10
        Ftv = 0.0;
    else
        Ftv = 8.0 * (rv0^2) * s2 * (cs + rt) / denom_v;
    end
    
    denom_h = cs * rt;
    if abs(denom_h) < 1e-10
        Fth = 0.0;
    else
        Fth = -8.0 * (rh0^2) * s2 * (cs + rt) / denom_h;
    end
    
    % Shadowing terms at nadir
    denom_st0v = abs(1.0 + 8.0 * rv0 / (cs * Ftv))^2;
    denom_st0h = abs(1.0 + 8.0 * rv0 / (cs * Fth))^2;
    
    if abs(denom_st0v) < 1e-10
        St0v = 1.0;
    else
        St0v = 1.0 / denom_st0v;
    end
    
    if abs(denom_st0h) < 1e-10
        St0h = 1.0;
    else
        St0h = 1.0 / denom_st0h;
    end
    
    % Compute series sums
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    temp1 = 1.0;
    
    ks_cs = ks * cs;
    ks_cs_sq = ks_cs^2;
    exp_ks_cs_sq = exp(-ks_cs_sq);
    
    for n = 1:min(n_terms, length(spectra))
        fn = double(n);
        temp1 = temp1 * (1.0 / fn);
        
        a0 = ks_cs^(2.0 * fn);
        weight = spectra(n);
        
        sum1 = sum1 + temp1 * a0 * weight;
        
        % Term for vertical polarization
        term_v = abs(Ftv / 2.0 + 2.0^(fn + 2.0) * rv0 / cs * exp_ks_cs_sq)^2;
        sum2 = sum2 + temp1 * a0 * term_v * weight;
        
        % BUG FIX 2: Term for horizontal polarization (use rh0, not rv0)
        term_h = abs(Fth / 2.0 + 2.0^(fn + 2.0) * rh0 / cs * exp_ks_cs_sq)^2;
        sum3 = sum3 + temp1 * a0 * term_h * weight;
    end
    
    % Shadowing terms
    if abs(sum2) < 1e-10
        Stv = 0.0;
    else
        Stv = 0.25 * abs(Ftv)^2 * sum1 / sum2;
    end
    
    if abs(sum3) < 1e-10
        Sth = 0.0;
    else
        Sth = 0.25 * abs(Fth)^2 * sum1 / sum3;
    end
    
    % Transition functions
    if abs(St0v) < 1e-10
        Tfv = 0.0;
    else
        Tfv = 1.0 - Stv / St0v;
    end
    
    if abs(St0h) < 1e-10
        Tfh = 0.0;
    else
        Tfh = 1.0 - Sth / St0h;
    end
    
    % Ensure non-negative
    Tfv = max(0.0, real(Tfv));
    Tfh = max(0.0, real(Tfh));
end
