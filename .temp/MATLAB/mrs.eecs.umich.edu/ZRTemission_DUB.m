%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 12.6: ZRT Emission Model for a layer with distinct upper boundary

%Description: Code computes emissivity of a snow layer above ground. See
%Fig. 12-17. The top boundary is assumed to be planar, but the lower
%boundary can be rough.

%Input Variables: 
    %eps2: dielectric constant of middle layer
    %eps3: dielectric constant of bottom medium
    %theta_i: incidence angle in air (deg)
    %f: frequency (GHz)
    %s: rms height of lower boundary (m)
    %l: correlation length of lower boundary (m)
    %a: Single scattering albedo 
    %d: middle layer thickness (m)
    %kappa_e: extinction coefficient of middle layer (Np/m)
    
%Output Products:
    %e_v(theta_i): v-polarized emissivity
    %e_h(theta_i): h-polarized emissivity
    
%Book Reference: Section 12-6

%Matlab code

function [e_v e_h] = ZRTemission_DUB(eps2,eps3, theta_i,f, s, l, a, d, kappa_e) 

eps1 = 1;
theta1 = theta_i*pi/180; %transform to radians
n2 = sqrt(eps2); %calc index of refraction
n1 = 1.0;          

theta2 = asin(abs(n1/n2) * sin(theta1)); % incidence angle in medium 2 (refracted)

%-- calculate reflectivies at the two interfaces

%- 12 interface
[rhoh rhov gammah12 gammav12 tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta_i);

%- 23 interface - use I2EM model 
[ev23 eh23] = Calc_emissivity_I2EMmodel(f, s, l, theta2*180/pi, eps3/eps2, 1);
gammav23 = 1- ev23;
gammah23 = 1- eh23;


trans = exp(-kappa_e * d /cos(theta2)); %extinction coefficient inside medium.

e_v = ( (1 - gammav12)/(1 - gammav12*gammav23*trans^2))* ...
    ( (1+gammav23*trans)*(1-a)*(1 - trans)+ (1 - gammav23)*trans);

e_h = ( (1 - gammah12)/(1 - gammah12*gammah23*trans^2))* ...
    ( (1+gammah23*trans)*(1-a)*(1 - trans)+ (1 - gammah23)*trans);

end
