%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.6: PRISM-1 Inverse Model

%Description: Given sigma_0_vv, sigma_0_hh, and sigma_0_hv, the inverse
%model computes eps, s, and mv (assuming loamy soil)

%Input Variables:
    %sigma_0_vv (dB), sigma_0_hh (dB), and sigma_0_hv (dB)
    %theta: Incidence angle (degrees)
    %f: frequency (GHz)
    
%Output Product:
    %eps: absolute value of complex dielectric constant of soil medium
    %s: rms height (m)
    %mv: volumetric moisture content (g/cm3)
    
%Book Reference: Section 10-5

%Matlab Code:

function [eps s mv] = PRISM1_InverseModel(sigma_0_vv,sigma_0_hh, sigma_0_hv, theta, f)

k = 2*pi*f ./0.3; %calculate the wave number 

theta_rad = theta*pi/180; %represent angle in radians

sigma_0_vv = 10.^(sigma_0_vv/10); %represent data in linear scale
sigma_0_hh = 10.^(sigma_0_hh/10);
sigma_0_hv = 10.^(sigma_0_hv/10);

p = sigma_0_hh ./ sigma_0_vv; %calculate the p-ratio
q = sigma_0_hv ./ sigma_0_vv; %calculate the q-ratio

gamma0 = 0:0.0001:1; % set gamma0 range of values (fine increments)

a1 = 2*theta_rad / pi;

%-- calculate the error function for different values of gamma0 and find
%the minimum

err = 1- a1.^(1./3./gamma0) .*(1- q./0.23./sqrt(gamma0)) - sqrt(p); 

abs_err = abs(err);
min_err = min(abs_err); %find the value of minimum error
gamma0_min = gamma0(abs_err == min_err);

eps = ( (1 + sqrt(gamma0_min)) ./(1- sqrt(gamma0_min)) ).^2;

s = -1/k .* log(1 - q ./(0.23*sqrt(gamma0_min)));


%-----------------------------------------------------------
%------ estimate mv based on dielectric model of soil (loamy)

t= 23; % temperature in degree C
S = 30.6 / 100; % fraction of Sand
C= 13.5 /100; % fraction of Clay

rho_b = 1.7; % soil bulk density g/cm3
f_hz = f * 1.0e9; % transform from GHz to Hz

beta1 = 1.27 - 0.519 * S - 0.152* C;  %eq: 4.68b
%beta2 = 2.06 - 0.928 * S - 0.255 * C; %eq: 44.68c 
alpha = 0.65; % eq: 4.68a

%Dielectric Constant of Pure Water
    
ew_inf = 4.9; % eq: E.15
ew_0 = 88.045 - 0.4147 * t + 6.295e-4 * t^2 + 1.075e-5 * t^3; %  
tau_w = (1.1109e-10 - 3.824e-12*t +6.938e-14*t^2 - 5.096e-16*t^3)/2/pi; %

epsrW = ew_inf +(ew_0-ew_inf)./(1 + (2*pi*f_hz*tau_w).^2); %


  %-- extract the value now---------------
  
mvv = 0:0.001:0.5;  
%-calculate error function
err2 = abs(eps.^(alpha) - 1- 0.66*rho_b - mvv.^beta1 * epsrW.^alpha + mvv);

min_err2 = min(err2); %find minimum error
mv = mvv(err2 == min_err2); %find mv corresponding to minimum error
end

