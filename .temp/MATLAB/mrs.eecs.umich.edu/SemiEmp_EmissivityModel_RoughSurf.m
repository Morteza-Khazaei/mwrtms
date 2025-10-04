%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code12.3: Semiempirical L-Band Emissivity Model for Rough Surface

%Description: Code computes polarized emissivity for a random surface using
%the parametric model equations developed by Lawrence et al. (2013). Model
%is applicable to L-band only.

%Input Variables:
    %eps: dielectric constant of emitting medium
    %s: rms heigh (m)
    %l: correlation length (m)
    %theta: Incidence angle (degrees)
   
%Output Product: 
    %e_v and e_h: emissivity for v and h polarizations.
    
%Book Reference: Section 12-3.2 and equations 12.11, 12.13-12.14

%Matlab Code:

function [e_v e_h] = SemiEmp_EmissivityModel_RoughSurf(eps, s,l,theta)

s = s * 100; % transform to cm
l = l * 100;

thetar = theta * pi/180;

Zs = s^2 /l; %cm

if Zs <= 1.235,
    h_pr = 2.615*(1 - exp(-Zs/2.473));
else
    h_pr = 1.0279;
end
Q = 0.1771 * h_pr;

n_v = 1.615 *(1-exp(-h_pr/0.359)) - 0.238; %for v-pol
n_h = 0.767* h_pr - 0.099; % for h-pol 


%--call function to calculate reflectivity of planar boundary!
[rhoh rhov gammah gammav tauh tauv Th ... 
        Tv] = ReflTransm_PlanarBoundary(1, eps, theta);

gamma_rs_v = ( (1-Q)*gammav + Q* gammah) * exp(-h_pr * (cos(thetar))^(n_v));

gamma_rs_h = ( (1-Q)*gammah + Q* gammav) * exp(-h_pr * (cos(thetar))^(n_h));

e_v = 1- gamma_rs_v;
e_h = 1- gamma_rs_h;


end


