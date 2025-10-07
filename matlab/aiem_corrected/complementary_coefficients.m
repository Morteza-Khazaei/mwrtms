function F = complementary_coefficients(u, v, q, qslp, qfix, R, eps_r, si, sis, cs, css, sfs, csfs, pol, is_substrate)
% COMPLEMENTARY_COEFFICIENTS Compute complementary field coefficients
%
% Computes complementary field coefficients for AIEM multiple scattering term.
%
% Inputs:
%   u, v         - Integration variables
%   q            - Wave vector component
%   qslp         - Slope-dependent wave vector
%   qfix         - Fixed wave vector component
%   R            - Reflection coefficient (Rv, Rh, or Rhv)
%   eps_r        - Complex relative permittivity
%   si, sis      - sin(theta_i), sin(theta_s)
%   cs, css      - cos(theta_i), cos(theta_s)
%   sfs, csfs    - sin(phi_s), cos(phi_s)
%   pol          - Polarization: 'vv', 'hh', 'hv', 'vh'
%   is_substrate - true for substrate-side (fb*), false for air-side (fa*)
%
% Output:
%   F - Complementary field coefficient
%
% Bug fixes applied:
%   1. Use abs(css - qslp) instead of abs(real(css - qslp))

    % Wave vector components
    kxu = si + u;
    ksxu = sis * csfs + u;
    kyv = v;
    ksyv = sis * sfs + v;
    
    % BUG FIX: Use complex magnitude, not real part magnitude
    if abs(css - qslp) < 1e-10
        zx = 0.0;
        zy = 0.0;
    else
        zx = -ksxu / (css - qslp);
        zy = -ksyv / (css - qslp);
    end
    
    if abs(cs + qslp) < 1e-10
        zxp = 0.0;
        zyp = 0.0;
    else
        zxp = kxu / (cs + qslp);
        zyp = kyv / (cs + qslp);
    end
    
    % Compute field components based on polarization
    switch lower(pol)
        case {'vv', 'hh'}
            % Co-polarization components
            c1 = -csfs * (-1.0 - zx * zxp) + sfs * zxp * zy;
            c2 = -csfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - ...
                          cs * v * zyp - si * v * zx * zyp) + ...
                  sfs * (cs * u * zy + si * u * zxp * zy + q * si * zyp - ...
                         cs * u * zyp + si * v * zy * zyp);
            c3 = -csfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) + ...
                  sfs * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy);
            c4 = -css * sfs * (-si * zyp + cs * zx * zyp) - ...
                  csfs * css * (-cs - si * zxp - cs * zy * zyp) + ...
                  sis * (-cs * zx - si * zx * zxp - si * zy * zyp);
            c5 = -css * sfs * (-v * zx + v * zxp) - ...
                  csfs * css * (q + u * zxp + v * zy) + ...
                  sis * (q * zx + u * zx * zxp + v * zxp * zy);
            c6 = -css * sfs * (-u * zyp + q * zx * zyp) - ...
                  csfs * css * (v * zyp - q * zy * zyp) + ...
                  sis * (v * zx * zyp - u * zy * zyp);
            
            % Reflection coefficient factors
            rp = 1.0 + R;
            rm = 1.0 - R;
            a = rp / qfix;
            b = rm / qfix;
            
            if strcmpi(pol, 'vv')
                if ~is_substrate
                    % Air-side VV
                    F = b * (-rp * c1 + rm * c2 + rp * c3) + a * (rm * c4 + rp * c5 + rm * c6);
                else
                    % Substrate-side VV
                    F = a * (rp * c1 - rm * c2 - rp * c3 / eps_r) - ...
                        b * (rm * c4 * eps_r + rp * c5 + rm * c6);
                end
            else  % HH
                if ~is_substrate
                    % Air-side HH
                    F = -b * (-rp * c1 + rm * c2 + rp * c3) - a * (rm * c4 + rp * c5 + rm * c6);
                else
                    % Substrate-side HH
                    F = a * (-rp * c1 * eps_r + rm * c2 + rp * c3) + ...
                        b * (rm * c4 + rp * c5 + rm * c6 / eps_r);
                end
            end
            
        case {'hv', 'vh'}
            % Cross-polarization components
            b1 = -css * sfs * (-1.0 - zx * zxp) - sis * zy - csfs * css * zxp * zy;
            b2 = -css * sfs * (-cs * q - cs * u * zx - q * si * zxp - si * u * zx * zxp - ...
                               cs * v * zyp - si * v * zx * zyp) + ...
                  sis * (-cs * q * zy - q * si * zxp * zy + q * si * zx * zyp - ...
                         cs * u * zx * zyp - cs * v * zy * zyp) - ...
                  csfs * css * (cs * u * zy + si * u * zxp * zy + q * si * zyp - ...
                                cs * u * zyp + si * v * zy * zyp);
            b3 = -css * sfs * (si * u - q * si * zx - cs * u * zxp + cs * q * zx * zxp) - ...
                  csfs * css * (-si * v + cs * v * zxp + q * si * zy - cs * q * zxp * zy) + ...
                  sis * (-si * v * zx + cs * v * zx * zxp + si * u * zy - cs * u * zxp * zy);
            b4 = -csfs * (-si * zyp + cs * zx * zyp) + sfs * (-cs - si * zxp - cs * zy * zyp);
            b5 = -csfs * (-v * zx + v * zxp) + sfs * (q + u * zxp + v * zy);
            b6 = -csfs * (-u * zyp + q * zx * zyp) + sfs * (v * zyp - q * zy * zyp);
            
            % Reflection coefficient factors
            rp = 1.0 + R;
            rm = 1.0 - R;
            a = rp / qfix;
            b = rm / qfix;
            
            if strcmpi(pol, 'hv')
                if ~is_substrate
                    % Air-side HV
                    F = b * (rp * b1 - rm * b2 - rp * b3) + a * (rm * b4 + rp * b5 + rm * b6);
                else
                    % Substrate-side HV
                    F = a * (-rp * b1 + rm * b2 + rp * b3 / eps_r) - ...
                        b * (rm * b4 * eps_r + rp * b5 + rm * b6);
                end
            else  % VH
                if ~is_substrate
                    % Air-side VH
                    F = b * (rp * b4 + rm * b5 + rp * b6) - a * (-rm * b1 + rp * b2 + rm * b3);
                else
                    % Substrate-side VH
                    F = -a * (rp * b4 + rm * b5 + rp * b6 / eps_r) + ...
                         b * (-rm * b1 * eps_r + rp * b2 + rm * b3);
                end
            end
    end
end
