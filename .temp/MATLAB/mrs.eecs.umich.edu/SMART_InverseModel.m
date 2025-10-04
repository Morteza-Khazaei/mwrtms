%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
%Code 10.8: SMART  - Inverse Model

%Description: Given sigma_0_vv and sigma_0_hh, the model computes eps_pr
%and s.  

%Input Variables:
    %sigma_0_vv (dB)
    %sigma_0_hh (dB)
    %theta: Incidence angle (degrees)
    %f: Frequency (GHz)
    
%Output Products:
    %eps_pr: real part of the dielectric constant of soil medium (unitless)
    %s: rms height of surface roughness (m)
    
    
%Book Reference: Section 10-6

%Matlab Code: 

function [eps_pr s] = SMART_InverseModel(sig_0_vv, sig_0_hh,theta,f)


lmda = 30/f; % wavelength in cm

k = (2*pi / lmda); % calculate wavenumber in cm

theta = theta .* (pi/180.0); % incidence angle in radian

eps_pr = 1/3.36 /tan(theta) .*(14 *sig_0_vv - 11*sig_0_hh + 26.5 ... 
    -255*log10(cos(theta)) -130*log10(sin(theta)) -21*log10(lmda));

logks = -0.08333 *sig_0_vv + 0.1369 *sig_0_hh + 1.806495 + 0.4465*log10(cos(theta)) ...
    +3.3451 * log10(sin(theta)) - 0.375 * log10(lmda);

s = 10.^(logks)/k;


end

