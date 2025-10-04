function [svv,shh,shv] = rough_bks_OH(eps1,eps2,tai,ks)
% calculate surface backscattering following the empirical model of Oh et
% al. 1992 of bare soil
% ks = k*s, k is in eps1, s is rms height
% Validated range:
%   0.1 < ks < 6, 2.5 < kl < 20, 0.09 < mv < 0.31, 
%   incidence angle 20~70 degree for smooth surface, and 0~70 for rough
%   surface.
% Reference: 
%   Oh, Sarabandi and Ulaby, TGRS 30(2): 370-381
% 

Gm0 = Fresnel(eps1,eps2,0);
[Gmv,Gmh] = Fresnel(eps1,eps2,tai);
q = 0.23*sqrt(Gm0)*(1 - exp(-ks));
pr = 1 - (2*tai/pi)^(1/3/Gm0)*exp(-ks);
p = pr^2;
g = 0.7*(1 - exp(-0.65*(ks)^1.8));

svv = g*cos(tai)^3/pr*(Gmv + Gmh);
shh = p*svv;
shv = q*svv;

end