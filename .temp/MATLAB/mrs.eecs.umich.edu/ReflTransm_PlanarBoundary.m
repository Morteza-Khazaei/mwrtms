%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 2.3: Oblique Reflection and Transmission @ Planar Boundry
%Description: Code computes the reflection coefficients, transmission
%coefficients, reflectivities and transmissivities for incidence in
%medium (medium 1) upon the planar boundary of a lossless or lossy
%medium (medium 2) at any incidence angle, for both h and v polarizations

%Input Variables:
    %eps1: eps1r -j*eps1i: relative dielectric constant of medium 1
    %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
    %theta1d: incidence angle in medium 1 in degrees

%Output Products:
    %rhoh: reflection coefficient for h pol
    %rhov: reflection coefficient for v pol
    %gammah:reflectivity for h pol
    %gammav: reflectivity for v pol
    %tauh: transmission coefficient for h pol
    %tauv: transmission coefficient for v pol
    %Th: transmissivity for h pol
    %Tv:transmissivity for v pol
%Book Reference: Sections 2-7 & 2-8

%Example call: [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)
%Computes 

%MATLAB Code

function [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)
    
    theta1 = degtorad(theta1d);
    
    sin_theta2 = sqrt(eps1)/sqrt(eps2).*sin(theta1);
    cos_theta2 = sqrt(1 - sin_theta2.^2);
    
    rhoh = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos_theta2) ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos_theta2);
    rhov = (sqrt(eps1).*cos_theta2-sqrt(eps2).*cos(theta1)) ./ (sqrt(eps1).*cos_theta2 + sqrt(eps2).*cos(theta1));

    tauh = 1 + rhoh;
    tauv = (1 + rhov).*(cos(theta1)./cos_theta2);
       

    gammah = abs(rhoh).^2;
    gammav = abs(rhov).^2;
        
    Th = 1-gammah;
    Tv = 1-gammav;
    
  end