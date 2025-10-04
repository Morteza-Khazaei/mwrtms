%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.8: Optical Thickness - Cloudy Atmosphere

%Description: Code computes the optical thickness under cloudy
    %conditions for an atmosphereic layer between heights z1 and z2,
    %containing a cloud between z3 and z4 for propagation along
    %direction theta relative to the zenith (0<theta<80) at any frquency
    %0<f<1000 GHz, using a Standard Atmospheric Model.  It is assumed that
    %the cloud layer is partially or totally within the atmospheric layer
    %of interest (z1 to z2).

%Input Variables
    %T0: Sea-surface temperatre in  degrees C
    %P0: sea-surface total barometric pressure in mbar
    %rhov0: sea-surface water vapor density in g/m^3
    %f: frequency in GHz
    %z1: base height of atmospheric layer in km
    %z2: top height of atmospheric layer in km
    %z3: base height of the cloud layer in km
    %z4: top height of the cloud layer in km
    %theta: zenith angle in degrees
    %cloud type: "water" or "ice"
    %mv: cloud waer content in g/m^3
    
%Output Products
    %tau: optical thickness in Np
    %Tr: Transmissivity (unitless)
%Book Reference: Sections 7-7

%MATLAB Code:

function [tau Tr]=OpticalThickn_CloudyAtmos(T0, P0, rhov0, f, z1, z2, z3, z4, theta, type, mv)  

%-- calculate the optical thickness of clear atmosphere 

  [tau_g Tr] = OpticalThickn_ClearAtmos(T0,P0,rhov0,f,z1,z2,theta); % 

%-- calculate the optical thickness of the cloud region
  
    int_tol = 1.0e-4; % set tolerance for the Quad integration function

    t = quad(@(z) calcCloud(z, T0,f, type, mv), z3, z4, int_tol);

    tau_c = t.*sec(deg2rad(theta));
    
    tau = tau_c + tau_g; % total optical thickness
    Tr = exp(-tau);
end

function [term]=calcCloud(z, T0, f,type,mv)

np = length(z);
term = zeros(np,1);
T0 = T0 + 273; % Transform to Kelvin

for ii=1:np
    %Calculate T
        if((z(ii)>=0)&&(z(ii)<11))
            T =T0-6.5*z(ii);
        elseif((z(ii)>=11)&&(z(ii)<20))
            T = T0-6.5*11;
        elseif((z(ii)>=20)&&(z(ii)<=32))
            T = (T0-6.5*11) +(z(ii)-20);
        end
        
        T = T - 273; %back to Celcius
        
    if (strcmp(type,'water'))    
        term(ii) = ExtincCoef_WaterCloud(T,mv,f); % in dB/km
    elseif (strcmp(type,'ice'))
        term(ii) = ExtincCoef_IceCloud(T,mv,f); % in dB/km
    else
        error('Invalid Cloud Type')
    end
    term(ii) = term(ii) / 4.34; % Transform to Np
    
end
end