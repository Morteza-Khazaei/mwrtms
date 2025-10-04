%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 12.7: Coherent Emissivity of a two-layer structure

%Description: Code computes emissivity of a perfectly homogeneous layer
%(such as ice or oil) above another homogeneous medium, bounded by
%perfectly smooth planar surfaces.

%Input Variables: 
    %eps2: dielectric constant of middle layer
    %eps3: dielectric constant of bottom layer
    %theta_i: incidence angle in air (deg)
    %d: layer thickness (m)
    %f: frequency (GHz)
    
%Output Products:
    %e_v_coh(theta_i): v-polarized emissivity
    %e_h_coh(theta_i): h-polarized emissivity
    
%Book Reference: Section 12-12.1

%Matlab code

function [e_v_coh e_h_coh] = CohEmissivity_TwoLayer(eps2, eps3, theta_i, d, f) 

eps1 = 1;
theta1 = theta_i*pi/180; %transform to radians
n2 = sqrt(eps2); %calc index of refraction
n1 = 1.;

theta2 = asin(abs(n1/n2) * sin(theta1)); % incidence angle in medium 2 (refracted)

k = 2*pi *f/ 0.3; %wavenumber 
%-- calculate reflectivies at the two interfaces

%- 12 interface
[rhoh12 rhov12 gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta_i);

%- 23 interface
[rhoh23 rhov23 gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps2, eps3, theta2*180/pi);


alpha = k* imag(sqrt(eps2));
beta = k * real(sqrt(eps2));

gamma = 1i * alpha +  beta;

rho_v = (rhov12 + rhov23 * exp(-2*1i*gamma*d*cos(theta2))) ./ ...
    (1 + rhov12 * rhov23 * exp(-2*1i*gamma*d*cos(theta2))); 

rho_h = (rhoh12 + rhoh23 * exp(-2*1i*gamma*d*cos(theta2))) ./ ...
    (1 + rhoh12 * rhoh23 * exp(-2*1i*gamma*d*cos(theta2))); 

e_v_coh = 1- (abs(rho_v)^2);
e_h_coh = 1- (abs(rho_h)^2);

end
