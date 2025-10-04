%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

vi = 0:0.01:1; % volume fraction of inclusion

eps_h = 1 - 1i*0; % complex dielectric constant of host material
eps_i = 10 - 1i* 1; %complex dielectric constant of inclusion material

shape = 2; % spheres

np = length(vi);

for ii = 1: np
    [eps_m(ii)] = TVBmodel_HeterogeneousMix(eps_i, eps_h, shape, vi(ii));
end

figure(1)
subplot(1,2,1)
plot(vi, real(eps_m))
grid

subplot(1,2,2)
plot(vi, abs(imag(eps_m)))
grid
