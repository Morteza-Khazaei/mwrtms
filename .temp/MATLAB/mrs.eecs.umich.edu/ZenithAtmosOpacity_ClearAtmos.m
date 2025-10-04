%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.9: Zenith Atmospheric Opacity - Clear Atmosphere

%Description: Using a Standard Atmosphere model, code computes opacity of
    %the entire atmosphere under clear-sky conditions, along the zenith
    %from sea-surface height z1=0 to z2=32km, at any
    %frequency 0<f<1000 GHz.

%Input Variables
    %T0: sea-surface temperatre in degree C
    %P0: sea-surface total barometric pressure in mbar
    %rhov0: sea-surface water vapor density in g/m^3
    %f: frequency in GHz
    
%Output Products:
    %tau0: zenith opacity in Np
    
%Book References: Sections 8-3

%MATLAB Code:


function [tau0]= ZenithAtmosOpacity_ClearAtmos(T0,P0,rhov0,f)

z1 = 0; % start height (sea level)
z2 = 32; % end height (max)
theta=0; % zenith incidence at 0 degrees

[tau0 Tr]= OpticalThickn_ClearAtmos(T0, P0, rhov0, f, z1, z2, theta); 
%- Tr is ignored here.

end