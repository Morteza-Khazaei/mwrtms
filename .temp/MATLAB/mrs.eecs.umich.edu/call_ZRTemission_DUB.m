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

f = 1.4; %freq (GHz)

s = 0.05; %rms height (m)
l = 0.2; %correlation length (m)

theta_i = 10; 

d = 0:0.01:2; %layer thickness in m

np = length(d); 

kappa_e = .5; %Np/m
a = 0.2; %single scattering albedo.

for n = 1:np
[e_v(n) e_h(n)] = ZRTemission_DUB(eps2,eps3, theta_i,f, s, l, a, d(n), kappa_e);
end

figure(1)
plot(d, e_v, '-r', d, e_h,'-b')
xlabel('thickness (m)')
ylabel('emissivity')
%axis([0 2 0.3 1])
grid

%

% % theta_i = 0:2:88; 
% % 
% % d = 0.2; %layer thickness in m
% % 
% % np = length(theta_i); 
% % 
% % kappa_e = .5; %Np/m
% % a = 0.2; %single scattering albedo.
% % 
% % for n = 1:np
% %     theta_i(n)
% % [e_v(n) e_h(n)] = ZRTemission_DUB(eps2,eps3, theta_i(n),f, s, l, a, d, kappa_e);
% % end
% % 
% % figure(1)
% % plot(theta_i, e_v, '-r', theta_i, e_h,'-b')
% % xlabel('Incidence Angle (deg)')
% % ylabel('emissivity')
% % %axis([0 2 0.3 1])
% % legend('e_v', 'e_h', 3)
% % grid

