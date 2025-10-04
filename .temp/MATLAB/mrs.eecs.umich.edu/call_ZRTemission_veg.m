%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

s = 0.05; %rms height (m)
l = 0.4; % correlation length (m)
f = 1.4; %freq in GHz
theta = 0:2:80;
np = length(theta);

eps = (9 - 1i*0.9)*2;
kappa_e = 0.2; % extinction rate Np/m
d = 2.0; %vegetation layer thickness (m)
a = 0.2; %albedo


for n = 1:np
    [e_v(n) e_h(n)] = ZRTemission_veg(eps, s,l,f,theta(n), a, kappa_e, d);
end

figure(1)

plot(theta, e_v, '-r', theta, e_h,'-b')
xlabel('\theta (deg)')
ylabel('Emissivity of Vegetation Layer')
legend('e_v', 'e_h', 3)
grid
