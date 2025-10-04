%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

rad = 0:0.001: 0.01; %particle radius in m
freq_GHz = 10.0 ; % freq in GHz

%-- permittivities of background and sphere

eps_bk = 1;  % real part of permittivity of background
eps_b =  eps_bk + 1i .* 0.0; % keep the background medium lossless

eps_sp = 3.1 + 1i * 2.1; % complex permittivity of spherical particle

%--------------

np = length(rad);

Es = zeros(np,1); Ea= Es; Eb=Es; Ee=Es;

eta_s_r = zeros(np,1); eta_a_r = eta_s_r; eta_e_r=eta_s_r; eta_b_r=eta_s_r;

for ii = 1:np
 [Es(ii) Ea(ii) Ee(ii) Eb(ii) eta_s_r(ii) eta_a_r(ii) eta_e_r(ii) eta_b_r(ii)] = ...
    Mie_Rayleigh_ScatteringOfSpheres(rad(ii), freq_GHz, conj(eps_sp), conj(eps_b)); 
end


figure(1)

subplot(2,2,1)
plot(rad, Es, rad, eta_s_r)
xlabel('radius (m)')
ylabel('Scattering Efficiency')
legend('Mie', 'Rayleigh')
grid

subplot(2,2,2)
plot(rad, Ea, rad, eta_a_r)
xlabel('radius (m)')
ylabel('Absorption Efficiency')
legend('Mie', 'Rayleigh')
grid

subplot(2,2,3)
plot(rad, Ee, rad, eta_e_r)
xlabel('radius (m)')
ylabel('Extinction Efficiency')
legend('Mie', 'Rayleigh')
grid
subplot(2,2,4)
plot(rad, Eb, rad, eta_b_r)
xlabel('radius (m)')
ylabel('Backscattering Efficiency')
legend('Mie', 'Rayleigh')
grid

