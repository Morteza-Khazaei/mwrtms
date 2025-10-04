%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.10: Zenith Atmospheric Opacity - Cloudy Atmosphere

%Description: Using a Standard Atmosphere model, code computes
    %opacity of the entire atmosphere when containing a cloud layer
    % extending between heights z3 and z4, along the zenith
    %from sea-surface height z1=0 to z2=32km, at any
    %frequency 0<f<50 GHz.

%Input Variables
    %T0: sea-surface temperatre in 0C
    %P0: sea-surface total barometric pressure in mbar
    %rhov0: sea-surface water vapor density in g/m^3
    %f: frequency in GHz
    %z3: base height of cloud layer in km
    %z4: top height of cloud layer in km
    %cloud type: water or ice
    %mv: cloud water content in g/m^3
    
%Output Products:
    %Tau0: zenith opacity in Np
    % Tr: Transmissivity
    
%Book References: Sections 8-3 and 8-7

%MATLAB Code:


function [tau0 Tr]=ZenithAtmosOpacity_CloudyAtmos(T0,P0,rhov0,f,z3,z4,type,mv)
z1=0;
z2=32;
theta=0;

[tau0 Tr]= OpticalThickn_CloudyAtmos(T0, P0, rhov0, f, z1, z2, z3, z4, theta, type, mv);

end