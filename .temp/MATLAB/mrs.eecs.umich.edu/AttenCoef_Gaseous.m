%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.3: Gaseous Attenuation Coefficient

%Description: Code computes the total attenuation (absorption) coefficient
    %due to Oxygen and water vapor at any frequency 1<f<1000GHz,
    %temperature -100<t<50C, pressure 10^-5 mbar < P < 1013 mbar,
    %and water vapor density 0<rho_0<20g/m^3.

%Input Variables
    %T: temperature in in degree C
    %P: total barometric pressure in mbar
    %rho_0: water vapor density in g/m^3
    %f: frequency in GHz
    
%Output Products
    %kappag: gaseous attenuation coefficient in dB/km
    
%Book Reference: Section 8-2.3

%MATLAB Code:

function [kappag]=AttenCoef_Gaseous(T,P,rho_0,f)

%-- calculate atten. due to water vapor
   [kappaH2O]=AttenCoef_WaterVapor(T,P,rho_0,f);

% calculate atten. due to oxygen    
   [kappaO2]=AttenCoef_Oxygen(T,P,rho_0,f);
   
% total atten. in dB/km   
   kappag=kappaH2O+kappaO2;
end