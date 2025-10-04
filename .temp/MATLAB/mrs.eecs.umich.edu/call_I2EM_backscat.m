%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
% This is a computationally demanding program, so it may take several
% minutes for it to execute. 


clear;

% fr = 1.5; % freq in GHz
% eps = 15.3 -1i* 3.7; % dielectric constant
fr = 4.75; 
eps = 15.2 -  1i* 2.12;

sig1 = 0.9 / 100 ; % rms height in m
L1 = 8.2 / 100; % correl length in m

sp = 1; xx= 1.5; % selection of correl func

theta = 10:1:80; %incidecne angle
np = length(theta);

sigma_0_hv = zeros(np,1);
sigma_0_vv = sigma_0_hv;
sigma_0_hh = sigma_0_hv;

for n = 1: np
    [sigma_0_vv(n) sigma_0_hh(n) sigma_0_hv(n)] = I2EM_Backscatter_model(...
    fr, sig1, L1, theta(n), eps, sp, xx);

end

figure(1)
plot(theta, sigma_0_vv,'-r', theta, sigma_0_hh, '-b', theta, sigma_0_hv, '-k')
xlabel('Incidence angle  \theta (deg)')
ylabel('\sigma^o (dB)')
legend('VV', 'HH', 'HV',4)
% axis([0 80 -40 10])
grid


