%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.7: SMART  - Forward Model

%Description: Code computes sigma_0_vv and sigma_0_hh, given surface 
%parameters.

%Input Variables:
    %eps = eps' - j eps'': Complex dielectric constant of the scattering
    %medium
    %theta: Incidence angle (degrees)
    %s: rms height (cm)
    %f: Frequency (GHz)
    
%Output Products:
    %sigma_0_vv (dB)
    %sigma_0_hh (dB)
    
%Book Reference: Section 10-6

%Matlab Code: 

function [sig_0_vv sig_0_hh] = SMART_ForwardModel(eps,theta,s,f)

eps_pr = real(eps); 

lmda = 30/f; % wavelength in cm

ks = s .* (2*pi / lmda); % calculate roughness parameter

theta = theta .* (pi/180.0); % incidence angle in radian

sig_0_hh = 10.^(-2.75) .* cos(theta).^(1.5) ./sin(theta).^5 .*lmda.^(0.7) ...
    .* (ks*sin(theta)).^(1.4) .* 10.^(0.028 * eps_pr .* tan(theta)); 

sig_0_vv = 10.^(-2.35) .* cos(theta).^(3) ./sin(theta).^3 .*lmda.^(0.7) ...
    .* (ks*sin(theta)).^(1.1) .* 10.^(0.046 * eps_pr .* tan(theta));

sig_0_vv = 10*log10(sig_0_vv);
sig_0_hh = 10*log10(sig_0_hh);


end

