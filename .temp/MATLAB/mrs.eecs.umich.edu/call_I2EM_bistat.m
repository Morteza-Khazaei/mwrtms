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

% sig = 0.246/100; % rms height (m)
% L = 0.02578; % correlation length (m)

% sig = 0.663/100; % rms height (m)
% L = 0.0506; % correlation length (m)

sig = 0.926/100; % rms height (m)
L = 0.0563; % correlation length (m)

sp = 1; % exponential correlation function
xx = 1.5; 

fr = 10; %frequency (GHz)

er = 3; % dielectric constant

thi = 20:.1:70; % incidence angle
phs = 0.0; 
np = length(thi);

sigma_0_vv = zeros(np,1);
sigma_0_hh = zeros(np,1);


for n = 1: np
 %-- using the I2EM bistatic code in specular
    ths = thi(n);
    [sigma_0_vv(n) sigma_0_hh(n)] = I2EM_Bistat_model(fr, sig, L, thi(n), ths,phs, er, sp,xx);
end


figure(1)
plot(thi, sigma_0_vv,'-r', thi, sigma_0_hh, '-b')
xlabel('Incidence angle (deg)')
ylabel('\sigma_0 (dB)')
legend('vv', 'hh',4)
axis([20 70 -50 10])
grid

%----------------------------------------------
%%
sig = 0.0955/100;
L = 0.477/100;
f = 10; 
thi = 45;
ths = 45;
er = 62;
sp = 2;

phs = 0:1:175;
np = length(phs);
sigma_0_vv = zeros(np,1);
sigma_0_hh = zeros(np,1);

for n = 1: np
 %-- using the I2EM bistatic code in specular
    [sigma_0_vv(n) sigma_0_hh(n)] = I2EM_Bistat_model(fr, sig, L, thi, ths,phs(n), er, sp,xx);
end


figure(2)
plot(phs, sigma_0_vv,'-b', phs, sigma_0_hh, '-r')
xlabel('Azimuth angle (deg)')
ylabel('\sigma_0 (dB)')
legend('vv', 'hh',4)
axis([0 175 -55 10])
grid






