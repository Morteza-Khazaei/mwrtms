%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.5: Ice Cloud Attenuation Coefficient

%Description: Code computes the extinction coefficient (scattering plus
    %absorption) for an ice cloud at any frequency 0 < f < 70 GHz,
    %temperature -50 < t < 0 C, and coud water content 0 < mv < 5g/m^3.

%Input Variables
    %T: temperature in degree C
    %mv: cloud water content in g/m^3
    %f: frequency in GHz
    
%Output Products
    %kappaic: cloud attenuation coefficient in dB/km
    
%Book Reference: Section 8-7.3

%MATLAB Code:

function [kappaic]= ExtincCoef_IceCloud(T,mv,f)

    lmda_0 = 30 ./f; % wavelength in cm

    %Dielecric Constant of pure ice
    [epsr epsi]=RelDielConst_PureIce(T,f);
    
    epsic = epsr - 1i*epsi;
    Nic2 = imag( -(epsic-1)./(epsic+2));

    kappaic=   0.434 * 6*pi/lmda_0 * mv * Nic2;
end