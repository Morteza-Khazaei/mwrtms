%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 6.2: Emissivity of Two-Layer Composite

%Description: Code computes the emissivity for emission by a two layer
    %lossy composite into a lossless medium (medium 1), for
    %both h and v polarizations.
    
%Input Variables:
    %eps1r: relative dielectric constant of medium 1
    %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
    %eps3 = eps3r-j*eps3i: relative dielectric constant of medium 3
    %d: vertical thickness of middle layer in meters
    %theta1: emission (transmission) angle in medium 1 in degrees
    %f: frequency in GHz
    
%Output Products:
    %eh: emissivity for h pol
    %ev: emissivity for v pol
%Book Reference: Section 6-7.5

%MATLAB Code

function [eh ev] = emissivity_TwoLayerComposite(eps1r, eps2, eps3, d, theta1, f)
    theta1 = degtorad(theta1);    
    eps1 = eps1r;
    
    lam2 = 1i*20*pi/3*f*sqrt(eps2);
    
    theta2 = acos((1-(sqrt(eps1)/sqrt(eps2).*sin(theta1)).^2).^(1/2));
    theta3 = acos((1-(sqrt(eps1)/sqrt(eps3).*sin(theta1)).^2).^(1/2));
    
    rho12h = (sqrt(eps1).*cos(theta1)-sqrt(eps2).*cos(theta2))...
            ./ (sqrt(eps1).*cos(theta1) + sqrt(eps2).*cos(theta2));
    rho12v = (sqrt(eps1).*cos(theta2)-sqrt(eps2).*cos(theta1))...
            ./ (sqrt(eps1).*cos(theta2) + sqrt(eps2)*cos(theta1));
    
    rho23h = (sqrt(eps2).*cos(theta2)-sqrt(eps3).*cos(theta3))...
            ./ (sqrt(eps2).*cos(theta2) + sqrt(eps3).*cos(theta3));
    rho23v = (sqrt(eps2).*cos(theta3)-sqrt(eps3).*cos(theta2))...
            ./ (sqrt(eps2).*cos(theta3) + sqrt(eps3)*cos(theta2));
    
    rhoh = (rho12h + rho23h.*exp(-2*lam2.*d.*cos(theta2)))...
            ./(1 + rho12h.*rho23h.*exp(-2*lam2.*d.*cos(theta2)));
    rhov = (rho12v + rho23v.*exp(-2*lam2.*d.*cos(theta2)))...
            ./(1 + rho12v.*rho23v.*exp(-2*lam2.*d.*cos(theta2)));

    gammah = abs(rhoh).^2;
    gammav = abs(rhov).^2;
    
    eh = 1-gammah;
    ev = 1-gammav;
    
  end