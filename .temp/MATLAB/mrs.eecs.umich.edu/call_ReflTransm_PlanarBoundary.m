%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps1 = 1 - 1i * 0.01; % relative dielectric constant of medium 1
eps2 = 50 -1i* 25; % relative dielectric constant of medium 2

theta1 = [0:1:90]; % incidence angle in medium 1

np = length(theta1);

for ii = 1:np
    [rhoh(ii) rhov(ii) gammah(ii) gammav(ii) tauh(ii) tauv(ii) Th(ii) Tv(ii)] = ReflTransm_PlanarBoundary(eps1, eps2, theta1(ii));
end

plot(theta1, abs(rhov), theta1, abs(rhoh))
xlabel('Incidence Angle (\theta_1)')
ylabel('|\rho|')
legend('|\rho_v|', '|\rho_h|',3)
grid