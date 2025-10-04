%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

epsr = [3.2:0.1:10]; % real part of relative dielectric constant
epsi = 0.5; % imaginary part

f = 10.0e9; % frequency in Hz

np = length(epsr);

for ii = 1: np
    [alpha(ii) beta(ii) mag_eta(ii) phase_eta(ii)]= SingleMedium_PropagConst(epsr(ii), epsi, f);
end

subplot(2,2,1)
plot(epsr, alpha)
grid

subplot(2,2,2)
plot(epsr, beta)
grid

subplot(2,2,3) 
plot(epsr, mag_eta)
grid

subplot(2,2,4)
plot(epsr, phase_eta)
grid


