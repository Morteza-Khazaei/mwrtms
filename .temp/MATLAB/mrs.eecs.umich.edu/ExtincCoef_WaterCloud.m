%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.4: Water Cloud Attenuation Coefficient

%Description: Code computes the extinction coefficient (scattering plus
    %absorption) for a water cloud at any frequency 0 < f < 50 GHz,
    %temperature 0 < t < 50 C, and cloud water content 0 < mv <5 g/m^3.

%Input Variables
    %T: temperature in degree C
    %mv: cloud water content in g/m^3
    %f: frequency in GHz
    
%Output Products
    %kappawc: cloud extinction coefficient in dB/km
    
%Book Reference: Section 8-7.3

%MATLAB Code:

function [kappawc]=ExtincCoef_WaterCloud(T,mv,f)

    lmda_0 = 30 ./f; % wavelength in cm
    
    %Dielecric Constant of pure water
    
    [epsr epsi] = RelDielConst_PureWater(T,f);
    
    epsw = epsr- 1i*epsi;
    
    k= imag(-(epsw-1)./(epsw+2));

    kappawc = 0.434 * 6*pi/lmda_0 * mv * k;
end