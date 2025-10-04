%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code12.1: Smooth Surface Emissivity

%Description: Code computes emissivity of a dielectric medium viewed from
%air(with eps1 = 1) at angle theta1 (Fig. 6-17)

%Input Variables:
    %eps: dielectric constant of emitting medium
    %theta1: Incidence angle in medium 1 (degrees)
   
%Output Product: 
    %e_v and e_h: emissivity for v and h polarizations.
    
%Book Reference: Section 6-7.2

%Matlab Code:

function [e_v e_h] = Emissivity_SmoothSurf(eps, theta1)

eps1r = 1; % dielectric constant of medium 1


[rhoh rhov gammah gammav tauh tauv Th ... 
        Tv] = ReflTransm_PlanarBoundary(eps1r, eps, theta1);

e_v = 1- gammav;
e_h = 1- gammah;
end


