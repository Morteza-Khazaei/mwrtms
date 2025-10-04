%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.6: Rain Attenuation Coefficient

%Description: Code computes the extinction coefficient (scattering plus
    %absorption) due to rain at any frequency 0<f<1000GHz,
    %rain rate 0<Rr<150mm/hr, and t>0C.

%Input Variables
    %Rr: rain rate in mm/hr
    %f: frequency in GHz
    %t: temperature in degree C
    
%Output Products
    %kapparain: rain attenuation coefficient in dB/km
    
%Book Reference: Sections 8-8.2

%MATLAB Code:

function [kapparain]=RainExtincCoef(Rr,f, t)

lmda = 0.3 ./ f; % wavelength in m

x_max = 2*pi ./ lmda *0.01; %maximum water droplet size 
x_min = 2*pi ./ lmda *1.0e-6;
[epsr epsi] = RelDielConst_PureWater(t,f); %calculate permittivity of water
eps = epsr -1i* epsi;

int_tol = 1.0e-2; % set tolerance for the Quad integration function

kapparain = quad(@(x) calcTermE(x, f, eps, Rr), x_min, x_max, int_tol);

kapparain = kapparain .* lmda.^3/8/ pi^2; % in Np/m

kapparain = kapparain * 4.3429*1000; % in dB/km

end

function [y] = calcTermE(x, f,eps, Rr)

np = length(x);

lmda = 0.3 ./f; %wavelength in meters
rad = x *lmda /(2*pi);

ext_mie= zeros(1, np); p_dist = ext_mie;

No = 8.0e6; 
b = 4100 * Rr^(-0.21);

for n = 1: np
   [Es Ea ext_mie(n) Eb eta_s_r eta_a_r eta_e_r eta_b_r] = ...
    Mie_Rayleigh_ScatteringOfSpheres(rad(n), f, eps, 1); 
    p_dist(n) = No .*exp(-b*2*rad(n)); %marchall palmer's function, normalized and scales!    
end

y = x.^2 .*p_dist .* ext_mie;
    
end
