%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 11.2: S2RT/R Backscattering from Rayleigh Layer with Distinct Upper
%Boundary

%Description: Code computes sigma_0_vv and sigma_0_hh for a weakly
%scattering Rayleigh layer with albedo a < 0.2. The layer is above a ground
%surface, and its upper and lower interfaces are characterized by the PRISM
%model (code 10.5).

%Input Variables:
    %eps2: complex dielectric constant of middle layer
    %eps3: complex dielectric constant of ground surface
    %f: frequency (GHz)
    %s1: rms height of top surface (m)
    %s3: rms height of ground surface (m)
    %a: Single-scattering albedo of Rayleigh layer (unitless)
    %kappa_e: extinction coefficient of the Rayleigh layer (Np/m)
    %d: thickness of the Rayleigh layer (m)
    %theta: Incidence angle (degrees)

%Output Variables: 
    %sigma_0_vv (dB)
    %sigma_0_hh (dB)
    
%Book reference: Section 11-2.2 and eq. 11.79 

%Matlab code:

function [sigma_0_vv sigma_0_hh] = S2RTR_DistinctUB(eps2,eps3,f,s1,s3,a,kappa_e,d, theta)

theta_r = theta *pi/180; %transform to radian


thetapr = asin(sqrt(1/real(eps2))*sin(theta_r)); % angle inside the layer #2
thetapr = thetapr *180/pi ;  % transform to degrees

costhetapr = sqrt( 1 - 1/real(eps2) .*sin(theta_r)^2);


kappa_s = a .* kappa_e ; %scattering coefficient

%---- call the PRISM-1 surface scattering model
%--  scattering due to top interface
[sig_s_vv sig_s_hh   tmp] = PRISM1_ForwardModel(eps2,theta,s1,f);
sig_12_vv  = 10.^(sig_s_vv/10); % transform to linear
sig_12_hh  = 10.^(sig_s_hh/10); 

%--  scattering due to bottom interface
    % note here we have assumed that PRISM1 model can be still used by
    % considering the ratio between the dielectric constants across the surface
    % boundary as the new dielectric constant input (i.e. used the dielectric
    % contrast as input).
[sig_s_vv sig_s_hh   tmp] = PRISM1_ForwardModel(eps3/eps2,thetapr,s3,f);
sig_23_vv  = 10.^(sig_s_vv/10); 
sig_23_hh  = 10.^(sig_s_hh/10); 


%-- calculate transmissivity in layer
tau = kappa_e .*d ./costhetapr;
T = exp(-tau);

%-- calculate reflectivity of upper interface
[t1 t2 t3 t4 t5 t6 Th_12 Tv_12] = ReflTransm_PlanarBoundary(1, eps2, theta);

%-- calculate reflectivity of lower interface
[t1 t2 gammah_23 gammav_23 t3 t4 t5 t6] = ReflTransm_PlanarBoundary(eps2, eps3, thetapr);

%-- calculate the total backscattering coefficient according to eq 11.79
sigma_0_vv = Tv_12.^2 .*( T.^2 .*sig_23_vv + 0.75*a * costhetapr .*(1- T.^2) ...
    .*(1 + gammav_23.^2 .* T.^2) + 3 * 2 * kappa_s * d * gammav_23 * T.^2) ...
    + sig_12_vv;

sigma_0_hh = Th_12.^2 .*( T.^2 .*sig_23_hh + 0.75*a * costhetapr .*(1- T.^2) ...
    .*(1 + gammah_23.^2 .* T.^2) + 3 * 2 * kappa_s * d * gammah_23 * T.^2) ...
    + sig_12_hh;



sigma_0_vv = 10* log10(sigma_0_vv);
sigma_0_hh = 10* log10(sigma_0_hh);


end