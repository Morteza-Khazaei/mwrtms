function [fvv, fhh, fhv, fvh] = kirchhoff_coefficients(Rv, Rh, k, theta_i, theta_s, phi_s, phi_i)
% KIRCHHOFF_COEFFICIENTS Compute Kirchhoff field coefficients
%
% Computes the Kirchhoff field coefficients for all polarizations
% using transitioned Fresnel reflection coefficients.
%
% Inputs:
%   Rv, Rh  - Transitioned reflection coefficients (V and H)
%   k       - Wavenumber (rad/m)
%   theta_i - Incident angle (radians)
%   theta_s - Scattered angle (radians)
%   phi_s   - Scattered azimuth (radians)
%   phi_i   - Incident azimuth (radians)
%
% Outputs:
%   fvv, fhh, fhv, fvh - Kirchhoff field coefficients

    % Trigonometric quantities
    si = sin(theta_i);
    cs = cos(theta_i);
    sis = sin(theta_s);
    css = cos(theta_s);
    sfs = sin(phi_s);
    csfs = cos(phi_s);
    
    % Cross-polarization coefficient
    Rhv = (Rv - Rh) / 2.0;
    
    % Surface slopes
    zxx = -(sis * csfs - si) / (css + cs);
    zyy = -(sis * sfs) / (css + cs);
    
    % Intermediate quantities
    d2 = sqrt((zxx * cs - si)^2 + zyy^2);
    
    % Field components
    hsnv = -(cs * csfs + si * (zxx * csfs + zyy * sfs));
    vsnh = css * csfs - zxx * sis;
    hsnh = -sfs;
    vsnv = zyy * cs * sis + css * (zyy * csfs * si - (cs + zxx * si) * sfs);
    
    hsnt = (-(cs^2 + si^2) * sfs * (-si + cs * zxx) + csfs * (cs + si * zxx) * zyy + si * sfs * zyy^2) / d2;
    hsnd = (-(cs + si * zxx) * (-csfs * si + cs * csfs * zxx + cs * sfs * zyy)) / d2;
    vsnt = ((cs^2 + si^2) * (-si + cs * zxx) * (csfs * css - sis * zxx) + ...
            css * sfs * (cs + si * zxx) * zyy - (csfs * css * si + cs * sis) * zyy^2) / d2;
    vsnd = -(cs + si * zxx) * (si * sis * zyy - css * (si * sfs - cs * sfs * zxx + cs * csfs * zyy)) / d2;
    
    % Kirchhoff field coefficients
    fhh = (1.0 - Rh) * hsnv + (1.0 + Rh) * vsnh - (hsnt + vsnd) * (Rh + Rv) * (zyy / d2);
    fvv = -((1.0 - Rv) * hsnv + (1.0 + Rv) * vsnh) + (hsnt + vsnd) * (Rh + Rv) * (zyy / d2);
    fhv = -(1.0 + Rv) * hsnh + (1.0 - Rv) * vsnv + (hsnd - vsnt) * (Rh + Rv) * (zyy / d2);
    fvh = -(1.0 + Rh) * hsnh + (1.0 - Rh) * vsnv + (hsnd - vsnt) * (Rh + Rv) * (zyy / d2);
end
