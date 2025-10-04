%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%---TSKY - Cloudy Atmosphere


function [TSKY] = TSKY_CloudyAtmos(T0,P0,rhov0,f,theta, z3,z4,type,mv)

    int_tol = 1.0e-2; % set tolerance for the Quad integration function
    T_extra = 2.7; % Temperature of Extra terrestrial radiation
    
   [tau0 Tr]=ZenithAtmosOpacity_CloudyAtmos(T0,P0,rhov0,f,z3,z4,type,mv);% calculate zenith opacity
    
   T_DN = quad(@(z) calcTerm(z, T0,P0,rhov0,f, theta,z3,z4,type,mv), 0, 32, int_tol); % Calc T_Downward
       
   T_DN = T_DN .*  sec(deg2rad(theta)); 
        
   TSKY = T_extra .* exp(-tau0 .* sec(deg2rad(theta))) + T_DN;
        
end

function [term]=calcTerm(z, T0,P0,rhov0,f, theta,z3,z4,type,mv)

    T0 = T0+273; % transform the temperature to Kelvin
        
    np = length(z); % size of vector z
    term = zeros(np,1);
    int_tol = 1.0e-4;    
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
    
    kappag = AttenCoef_Gaseous(T-273,P,rhov,f); % calc kappa_gaseous
    kappag = kappag /4.34; % transform the atten to Np
    
    if z(ii)<= z4 && z(ii)>=z3 
        if (strcmp(type,'water'))    
            kappac = ExtincCoef_WaterCloud(T-273,mv,f) / 4.34; 
        elseif (strcmp(type,'ice'))
            kappac = ExtincCoef_IceCloud(T-273,mv,f) / 4.34; 
        end    
    else
        kappac = 0;
    end
    kappa = kappag + kappac;
    
    %calculate optical thickness at height z(ii)
    
    %-- clear atmosphere:
    [tau_g Tr] = OpticalThickn_ClearAtmos(T0-273,P0,rhov0,f,0,z(ii),theta);
    
    %-- cloudy region
    if z(ii) > z4
        tau_c = quad(@(zz) calcCloud2(zz, T0-273,f, type, mv), z3, z4, int_tol); %calc kappa cloud    
    elseif z(ii) >=z3
          tau_c = quad(@(zz) calcCloud2(zz, T0-273,f, type, mv), z3, z(ii), int_tol); %calc kappa cloud    
        else 
            tau_c = 0;
    end
    tau_c = tau_c * sec(deg2rad(theta));
    
    tau_z = tau_g + tau_c;
    
    term(ii) = kappa * T * exp(-tau_z);  

    end
        
end

function [term]=calcCloud2(z, T0, f,type,mv)

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

