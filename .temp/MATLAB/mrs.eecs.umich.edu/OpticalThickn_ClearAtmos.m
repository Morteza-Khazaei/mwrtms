%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.7: Optical Thickness - Clear Atmosphere

%Description: Code computes the optical thickness under clear-sky
    %conditions for an atmosphereic layer between heights z1 and z2 for
    %propagation along direction theta relative to the zenith (0<theta<80)
    %at any ffrquency 0<f<1000 GHz, using a Standard Atmospheric Model.

%Input Variables
    %T0: Sea-surface temperature in degrees C
    %P0: sea-surface total barometric pressure in mbar
    %rhov0: sea-surface water vapor density in g/m^3
    %f: frequency in GHz
    %z1: base height of atmospheric layer in km
    %z2: top height of atmospheric layer in km
    %theta: zenith angle in degrees
    
%Output Products
    %tau: optical thickness in Np
    %Tr: Transmissivity (unitness)
%Book Reference: Sections 8-3

%MATLAB Code:

function [tau Tr]=OpticalThickn_ClearAtmos(T0, P0, rhov0, f, z1, z2, theta) 

    int_tol = 1.0e-4; % set tolerance for the Quad integration function

    t = quad(@(z) calcTerm(z, T0,P0,rhov0,f), z1, z2, int_tol);
    tau = t.*sec(deg2rad(theta)); % optical thickness in Np
    Tr = exp(-tau ); 
end
%
function [term]=calcTerm(z, T0,P0,rhov0,f)

    T0 = T0+273; % transform the temperature to Kelvin
        
    np = length(z); % size of vector z
    term = zeros(np,1);
        
    for ii = 1:np

    %Calculate T
        if((z(ii)>=0)&&(z(ii)<11))
            T=T0-6.5*z(ii);
        elseif((z(ii)>=11)&&(z(ii)<20))
            T = T0-6.5*11;
        elseif((z(ii)>=20)&&(z(ii)<=32))
            T = (T0-6.5*11) +(z(ii)-20);
        end
    
    %Calculate P
        P = P0*exp(-z(ii)/7.7);
    
    %Calculate rhov
        rhov = rhov0*exp(-z(ii)/2.1);
    
        T = T - 273; % transform back to C
        term(ii) = AttenCoef_Gaseous(T,P,rhov,f);
        term(ii) = term(ii) /4.34; % transform the atten to Np
        
%        z(ii), term(ii)
    end
        
end


