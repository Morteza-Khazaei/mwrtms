%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;


sig1 = 0.4/100; % rms height (m)
L1 = 0.07; % correlation length (m)

sig2 = 0.25/100; % rms height (m)
L2 = 0.03; % correlation length (m)

sig3 = 0.13/100; % rms height (m)
L3 = 0.015; % correlation length (m)
er = 9; % dielectric constant


% sig1 = 1./100; % rms height (m)
% L1 = 0.25; % correlation length (m)
% 
% sig2 = 0.7/100; % rms height (m)
% L2 = 0.11; % correlation length (m)
% 
% sig3 = 0.5/100; % rms height (m)
% L3 = 0.05; % correlation length (m)
% er = 5000000; % dielectric constant


fr = 4; %frequency (GHz)


thi = 0:1:60; % incidence angle
np = length(thi);

sigma_0_vv = zeros(np,1);
sigma_0_hh = zeros(np,1);


for n = 1: np
 %-- using the I2EM Multiscale code
    [sigma_0_vv(n) sigma_0_hh(n)] = Multiscale_I2EM_Backscatter(fr, ...
        sig1, L1,sig2, L2, sig3, L3, thi(n), er);
end


figure(1)
plot(thi, sigma_0_vv,'-r', thi, sigma_0_hh, '-b')
xlabel('Incidence angle (deg)')
ylabel('\sigma_0 (dB)')
legend('vv', 'hh',4)
%axis([10 80 -50 10])
grid

