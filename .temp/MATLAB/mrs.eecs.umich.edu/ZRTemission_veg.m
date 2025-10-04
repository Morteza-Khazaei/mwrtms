%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code12.5: ZRT emission model for vegetation

%Description: Code computes emissivity of a vegetation layer above a soil
%surface. 

%Input Variables:
    %eps: dielectric constant of emitting medium
    %s: rms heigh (m)
    %l: correlation length (m)
    %f: frequency (GHz)
    %theta: Incidence angle (degrees)
    %a: single scattering albedo
    %kappa_e: extinction coefficient of vegetation layer (Np/m)
    %d: depth of vegetation layer (m)
   
%Output Product: 
    %e_v and e_h: emissivity for v and h polarizations.
    
%Book Reference: Section 12-5.3

%Matlab Code:

function [e_v e_h] = ZRTemission_veg(eps, s,l,f,theta, a, kappa_e, d)

k = 2*pi*f/0.3; % wave number
thetar = theta * pi/180;

%--call function to calculate reflectivity of planar boundary!
[rhoh rhov gammah gammav tauh tauv Th ... 
        Tv] = ReflTransm_PlanarBoundary(1, eps, theta);


%--calculate coherent reflectivity
gamma_coh_v = gammav * exp(-(2*k*s* cos(thetar))^2);
gamma_coh_h = gammah * exp(-(2*k*s* cos(thetar))^2);


% emissivity of soil. Use I2EM emissivity code for rough surface
[e_v_soil e_h_soil] = Calc_emissivity_I2EMmodel(f, s, l, theta, eps, 1);

% attentuation through vegetation layer
veg_att = exp(-kappa_e*d * sec(thetar));

%total emissivity of the vegetation layer over rough soil
e_v = (1+ gamma_coh_v * veg_att) *(1-a)*(1- veg_att) + veg_att* e_v_soil;

e_h = (1+ gamma_coh_h * veg_att) *(1-a)*(1- veg_att) + veg_att* e_h_soil;


end


