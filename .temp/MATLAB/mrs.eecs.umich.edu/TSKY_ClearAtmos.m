%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 8.11: TSKY - Clear Atmosphere

%Description: Using a Standard Atmosphere model, code computes the
    %brightness temperature of the atmosphere as would be observed by a
    %ground-based upward-looking radiometer, under clear-sky conditions,
    %along direction theta relative to the zenith (0<theta<80degrees), from
    %sea-surface height z1=0 to z2=32km, at any frequency 1<f<1000 GHz.

%Input Variables:
    %T0: sea-surface temperature in degree C
    %P0: sea-surface total barometric presure in mbar
    %rhov0: sea-surface water vapor density in g/m^3
    %f: frequency in GHz
    %theta: zenith angle in degrees
    
%Output Products:
    %TSKY: Sky radiometric temperature in Kelvins
    
%Book Reference: Section 8-4

%MATLAB Code

function [TSKY] = TSKY_ClearAtmos(T0,P0,rhov0,f,theta)

    int_tol = 1.0e-2; % set tolerance for the Quad integration function
    T_extra = 2.7; % Temperature of Extra terrestrial radiation
    
    tau0 = ZenithAtmosOpacity_ClearAtmos(T0,P0,rhov0,f); % calculate zenith opacity

        T_DN = quad(@(z) calcTerm(z, T0,P0,rhov0,f, theta), 0, 32, int_tol); % Calc T_Downward
        T_DN = T_DN .*  sec(deg2rad(theta)); 
        TSKY = T_extra .* exp(-tau0 .* sec(deg2rad(theta))) + T_DN;
        
end

function [term]=calcTerm(z, T0,P0,rhov0,f, theta)

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
    
    kappa = AttenCoef_Gaseous(T-273,P,rhov,f); % calc kappa_gaseous
    kappa = kappa /4.34; % transform the atten to Np
    
    [tau_z Tr] = OpticalThickn_ClearAtmos(T0-273, P0, rhov0, f, 0, z(ii), theta);
    
    term(ii) = kappa * T * exp(-tau_z);  

    end
        
end

