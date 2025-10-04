%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 6.1: Emissivity of Homogeneous Medium with Planar Boundary

%Description: Code computes the emissivity for emission by a lossy medium
    %(medium 2) with a planar boundary into a lossless medium (medium 1), for
    %both h and v polarizations.
    
%Input Variables:
    %eps1r: relative dielectric constant of medium 1
    %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
    %theta1: emission (transmission) angle in medium 1 in degrees
    
%Output Products:
    %eh: emissivity for h pol
    %ev: emissivity for v pol
    
%Book Reference: Sections 6-7.2

%MATLAB Code

function [eh ev] = emissivity_PlanarBoundary(eps1r, eps2, theta1)
    theta1 = degtorad(theta1);
    eps1 = eps1r;
    
    theta2 = acos((1-(sqrt(eps1)/sqrt(eps2).*sin(theta1)).^2).^(1/2));
    rhoh = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos(theta2)) ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos(theta2));
    rhov = (sqrt(eps1).*cos(theta2)-sqrt(eps2).*cos(theta1)) ./ (sqrt(eps1).*cos(theta2) + sqrt(eps2).*cos(theta1));

    gammah = abs(rhoh).^2;
    gammav = abs(rhov).^2;
    
    eh = 1-gammah;
    ev = 1-gammav;
    
  end