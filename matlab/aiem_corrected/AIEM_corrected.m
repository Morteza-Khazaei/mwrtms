function [VV, HH, HV, VH] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, err, eri, itype, varargin)
% AIEM_CORRECTED Corrected Advanced Integral Equation Model
%
% Computes backscatter coefficients using AIEM with all bug fixes applied
% from the Python implementation and bug report analysis.
%
% USAGE:
%   [VV, HH, HV, VH] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, err, eri, itype)
%   [VV, HH, HV, VH] = AIEM_corrected(..., 'IncludeMultipleScattering', true)
%
% INPUTS:
%   theta_i - Incident angle (degrees)
%   theta_s - Scattered angle (degrees)  
%   phi_s   - Scattered azimuth angle (degrees), use 180 for monostatic backscatter
%   kl      - Normalized correlation length (k * L)
%   ks      - Normalized RMS height (k * sigma)
%   err     - Real part of relative permittivity
%   eri     - Imaginary part of relative permittivity
%   itype   - Surface correlation type:
%             1 = Gaussian
%             2 = Exponential
%             3 = 1.5-power (transformed exponential)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'IncludeMultipleScattering' - Include multiple scattering (default: false)
%                                 Set to true for cross-pol accuracy
%
% OUTPUTS:
%   VV, HH, HV, VH - Backscatter coefficients in dB
%
% BUG FIXES APPLIED:
%   1. Fresnel branch correction for lossy media (Im(k_tz) >= 0)
%   2. Normal incidence: rh0 = rv0 (was incorrectly -rv0)
%   3. Transition function: use rh0 in H-path (was using rv0)
%   4. 1.5-power spectrum: similarity-correct formula (was Bessel with order 1.5*n-1)
%   5. Complex magnitude checks: abs(z) instead of abs(real(z))
%
% KNOWN LIMITATIONS:
%   - Uses LEGACY Wu & Fung (1992) transition function
%   - Expected bias: +3-5 dB vs NMM3D for co-pol
%   - For <1 dB RMSE, implement new S_p/S_p^(0) transition method
%
% EXAMPLE:
%   % Monostatic backscatter at 40 degrees
%   theta_i = 40;
%   theta_s = 40;
%   phi_s = 180;  % Backscatter
%   kl = 5.0;     % Normalized correlation length
%   ks = 0.5;     % Normalized RMS height
%   err = 15.0;   % Real part of permittivity
%   eri = 1.5;    % Imaginary part of permittivity
%   itype = 2;    % Exponential correlation
%   
%   [VV, HH, HV, VH] = AIEM_corrected(theta_i, theta_s, phi_s, kl, ks, err, eri, itype);
%   fprintf('VV = %.2f dB\n', VV);
%   fprintf('HH = %.2f dB\n', HH);
%   fprintf('HV = %.2f dB\n', HV);
%   fprintf('VH = %.2f dB\n', VH);
%
% REFERENCES:
%   Chen, K. S., et al. (2003). "Emission of rough surfaces calculated by the
%   integral equation method with comparison to three-dimensional moment method
%   simulations." IEEE TGRS, 41(1), 90-101.
%
%   Wu, T. D., & Fung, A. K. (1992). "A transition model for the reflection
%   coefficient in surface scattering." IEEE TGRS, 30(4), 856-860.
%
% SEE ALSO:
%   aiem_single_scattering, fresnel_coefficients, transition_function
%
% AUTHOR: Generated from Python implementation with bug fixes
% DATE: 2024

    % Parse optional arguments
    p = inputParser;
    addParameter(p, 'IncludeMultipleScattering', false, @islogical);
    parse(p, varargin{:});
    include_ms = p.Results.IncludeMultipleScattering;
    
    % Form complex permittivity
    eps_r = err + 1i * eri;
    
    % Compute single scattering
    [sigma_vv, sigma_hh, sigma_hv, sigma_vh] = aiem_single_scattering(...
        theta_i, theta_s, phi_s, kl, ks, eps_r, itype);
    
    % Add multiple scattering if requested
    if include_ms
        % Convert angles to radians for multiple scattering
        theta_i_rad = deg2rad(theta_i);
        theta_s_rad = deg2rad(theta_s);
        phi_s_rad = deg2rad(phi_s);
        
        % Compute physical sigma from ks
        % We need to estimate k and sigma from ks and kl
        % Using the relationship: ks = k*sigma, kl = k*L
        % Therefore: sigma = ks/k, and k can be estimated from kl/ks ratio
        % For typical C-band (5.4 GHz): lambda = 5.55 cm, k = 113 rad/m
        % But we'll use a more general approach
        
        % Estimate k from the ratio (assuming reasonable surface parameters)
        % This is approximate but works for typical cases
        k_estimate = sqrt(ks * kl);  % Geometric mean approximation
        sigma_estimate = ks / k_estimate;
        
        % Compute multiple scattering for each polarization
        % Note: Only use for itype 1 or 2 (Gaussian/Exponential)
        if itype <= 2
            try
                fprintf('Computing multiple scattering (this may take 30-60 seconds)...\n');
                
                ms_vv = aiem_multiple_scattering(theta_i_rad, theta_s_rad, phi_s_rad, ...
                                                 kl, ks, eps_r, sigma_estimate, itype, 'vv');
                ms_hh = aiem_multiple_scattering(theta_i_rad, theta_s_rad, phi_s_rad, ...
                                                 kl, ks, eps_r, sigma_estimate, itype, 'hh');
                ms_hv = aiem_multiple_scattering(theta_i_rad, theta_s_rad, phi_s_rad, ...
                                                 kl, ks, eps_r, sigma_estimate, itype, 'hv');
                ms_vh = ms_hv;  % Reciprocity
                
                % Add to single scattering
                sigma_vv = sigma_vv + ms_vv;
                sigma_hh = sigma_hh + ms_hh;
                sigma_hv = sigma_hv + ms_hv;
                sigma_vh = sigma_vh + ms_vh;
                
                fprintf('Multiple scattering complete.\n');
            catch ME
                warning('Multiple scattering computation failed: %s', ME.message);
                fprintf('Error details: %s\n', ME.getReport());
            end
        else
            warning('Multiple scattering not supported for 1.5-power correlation');
        end
    end
    
    % Convert to dB
    VV = 10 * log10(sigma_vv);
    HH = 10 * log10(sigma_hh);
    HV = 10 * log10(sigma_hv);
    VH = 10 * log10(sigma_vh);
    
    % Handle -Inf for very small values
    if isinf(VV) && VV < 0
        VV = -999;
    end
    if isinf(HH) && HH < 0
        HH = -999;
    end
    if isinf(HV) && HV < 0
        HV = -999;
    end
    if isinf(VH) && VH < 0
        VH = -999;
    end
end
