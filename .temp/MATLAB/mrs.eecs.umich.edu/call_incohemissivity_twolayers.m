%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps2 = 3.2 -1i*0.001;
eps3 = 60- 1i*10;

theta_i = 10; 

d = 0:0.01:2; %layer thickness in m

np = length(d); 

kappa_e = .5; %Np/m
a = 0.2; %single scattering albedo.

for n = 1:np
[e_v_inc(n) e_h_inc(n)] = IncohEmissivity_TwoLayer(eps2,eps3, ...
    theta_i,a, d(n), kappa_e);
end

figure(1)
plot(d, e_v_inc, '-r', d, e_h_inc,'-b')
xlabel('thickness (m)')
ylabel('emissivity')
%axis([0 2 0.3 1])
grid