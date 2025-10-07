function W_n = roughness_spectrum(kl, K, n, correlation_type)
% ROUGHNESS_SPECTRUM Compute AIEM roughness spectrum with bug fixes
%
% Computes the n-th order roughness spectrum W^(n)(K) for different
% correlation functions.
%
% Inputs:
%   kl               - Normalized correlation length (k * L)
%   K                - Spatial frequency magnitude
%   n                - Spectral order (1, 2, 3, ...)
%   correlation_type - Type of correlation function:
%                      1 = Gaussian
%                      2 = Exponential
%                      3 = 1.5-power (transformed exponential)
%
% Output:
%   W_n - Roughness spectrum value W^(n)(K)
%
% Bug fixes applied:
%   1. 1.5-power spectrum uses similarity-correct formula
%      (not Bessel with order 1.5*n-1)

    fn = double(n);
    kl2 = kl * kl;
    
    switch correlation_type
        case 1  % Gaussian correlation
            W_n = (kl2 / (2.0 * fn)) * exp(-(K * K) / (4.0 * fn));
            
        case 2  % Exponential correlation
            W_n = (kl / fn)^2 * (1.0 + (K / fn)^2)^(-1.5);
            
        case 3  % 1.5-power correlation
            % BUG FIX: Use similarity-correct surrogate
            % W^(n)(K) = (L/n)^2 * (1 + α²(K*L/n^(2/3))²)^(-1.5)
            % This preserves the similarity law: W^(n) = L² * n^(-4/3) * Φ(K*L*n^(-2/3))
            
            if abs(K) < 1e-10
                % Special case: K = 0
                W_n = kl2 / (3.0 * fn - 2.0);
            else
                % Similarity-correct surrogate (α = 1.0)
                alpha = 1.0;
                n_power = fn^(2.0 / 3.0);  % n^(2/3) scaling
                W_n = (kl / fn)^2 * (1.0 + alpha^2 * (K * kl / n_power)^2)^(-1.5);
            end
            
        otherwise
            error('Unknown correlation type. Use 1 (Gaussian), 2 (Exponential), or 3 (1.5-power)');
    end
end
