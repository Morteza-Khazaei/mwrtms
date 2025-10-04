function [rv,rh] = rough_ref_WM(eps1,eps2,tai,ks)
% calculate reflectivity of rough interface using the semiempirical model
% by Wegmuller and Matzler, 1999
% ks = k*s, k is in eps1, s is rms height
% Validation range:
%   0.07 < ks < 27.5, 1~100GHz, 0~70 degree incidence, H and V pol
% Reference:
%   Wegmuller and Matzler, TGRS, 37(3): 1391-1395, 1999
% 

[rv0,rh0] = Fresnel(eps1,eps2,tai);
rh = rh0*exp(-ks^sqrt(0.10*cos(tai)));
degi = tai/pi*180;
if degi < 60
    rv = rh*cos(tai)^0.655;
else % 60~70
    rv = rh*(0.635 - 0.0014*(degi - 60));
end

end