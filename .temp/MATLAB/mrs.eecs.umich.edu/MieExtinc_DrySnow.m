%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 11.3: Mie Extinction of Dry Snow

%Description: Code computes the Mie absorption, scattering, and extinction
%coefficients of dry snow for specified snow density and ice particle
%radius.

%Input Variables:
    % rho_s: Snow density (g/cm3)
    % ri: Ice-particle radius (m)
    % f: frequency (GHz)
    % t: temperature (degree C)
    
%Output Products:
    % kappa_a: absorption coefficient (Np/m)
    % kappa_s: scattering coefficient (Np/m)
    % kappa_e: extinction coefficient (Np/m)
    % a: single-scattering albedo
    
%Book Reference: Section 11-15.1 and eq 11.112 and 11.113 with Q computed
%according to the Mie model of section 8-5.

%Matlab code

function [kappa_a kappa_s kappa_e a] = MieExtinc_DrySnow(rho_s,ri,f,t)

eps_b = 1; % dielectric constant of background
rho_i = 0.9167; %density of ice (g /cm3)

%- calculate relative dielectric constant of pure ice
[epsr epsi] = RelDielConst_PureIce(t,f);
eps_sp = epsr -1i * epsi;

%-- calculate the Mie efficiencies using Code 8.12
[Es Ea Ee Eb t1 t2 t3 t4] = ...
    Mie_Rayleigh_ScatteringOfSpheres(ri, f, eps_sp, eps_b); 

area = pi * ri^2;

Qs = Es * area;
Qa = Ea * area;
Qe = Ee * area;

Nv = rho_s ./rho_i ./(4/3 * pi .*ri^3); 

kappa_a = Nv .* Qa;
kappa_s = Nv .* Qs; 
kappa_e = Nv .* Qe;

a = kappa_s /kappa_e;

end

    