%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 10; %frequency (GHz)
s = 0.001; % rms surface height (m)

eps = 6.28 - 1i * 1.53; % corresponding to mv = 0.16

a = 0.15; %albedo must be < 0.2

kappa_e = 0.1 ; % extinction coefficient in Np/m

d = 10; %layer thickness in meters


theta = 5:1:80; %incidence angle (deg)
nt = length(theta);

%--calculate sigma_0 for the Rayleigh layer
for ii = 1:nt
[sig_0_vv(ii) sig_0_hh(ii) ] = S2RTR_DiffuseUB(eps,f,s,a,kappa_e,d, theta(ii));
end

figure(1)

plot(theta, sig_0_vv, theta, sig_0_hh)
xlabel('\theta (deg)')
ylabel('\sigma^0')
legend('vv', 'hh')
grid



