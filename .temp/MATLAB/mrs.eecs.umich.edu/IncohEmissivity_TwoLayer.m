%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 12.8: Incoherent Emissivity of a two-layer structure

%Description: Code computes incoherent emissivity of an inhomogeneous layer
%separated by air on top and a homogeneous medium on the bottom, with
%perfectly smooth parallel boundaries on both sides.

%Input Variables: 
    %eps2: dielectric constant of middle layer
    %eps3: dielectric constant of bottom layer
    %theta_i: incidence angle in air (deg)
    %a: Single scattering albedo
    %d: layer thickness (m)
    %kappa_e: extinction rate through layer (Np/m)
    
%Output Products:
    %e_v_inc(theta_i): v-polarized emissivity
    %e_h_inc(theta_i): h-polarized emissivity
    
%Book Reference: Section 12-12.2

%Matlab code

function [e_v_inc e_h_inc] = IncohEmissivity_TwoLayer(eps2,eps3, ...
    theta_i,a, d, kappa_e) 

eps1 = 1;
theta1 = theta_i*pi/180; %transform to radians
n2 = sqrt(eps2); %calc index of refraction
n1 = 1.0;

theta2 = asin(abs(n1/n2) * sin(theta1)); % incidence angle in medium 2 (refracted)

%-- calculate reflectivies at the two interfaces

%- 12 interface
[rhoh rhov gammah12 gammav12 tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta_i);

%- 23 interface
[rhoh rhov gammah23 gammav23 tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps2, eps3, theta2*180/pi);


trans = exp(-kappa_e * d /cos(theta2)); %extinction coefficient inside medium.

e_v_inc = ( (1 - gammav12)/(1 - gammav12*gammav23*trans^2))* ...
    ( (1+gammav23*trans)*(1-a)*(1 - trans)+ (1 - gammav23)*trans);

e_h_inc = ( (1 - gammah12)/(1 - gammah12*gammah23*trans^2))* ...
    ( (1+gammah23*trans)*(1-a)*(1 - trans)+ (1 - gammah23)*trans);

end
