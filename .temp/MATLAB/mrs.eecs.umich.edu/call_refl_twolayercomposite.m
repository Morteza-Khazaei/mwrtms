%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps2 = 3.2 - 1i * 0.2; % relative dielectric constant of medium 2
eps3 = 87 -1i* 58; % relative dielectric constant of medium 3

theta1 = [0:1:90]; % incidence angle in medium 1 (degrees)

d = 0.2 ; %layer thickness in m
f = 1.0; %frequency in GHz

np = length(theta1);

for ii = 1:np
    [rhoh(ii) rhov(ii) gammah(ii) gammav(ii)] = Refl_TwoLayerComposite(eps2, eps3, d, theta1(ii), f);
end

plot(theta1, abs(rhov), theta1, abs(rhoh))
xlabel('Incidence Angle (\theta_1)')
ylabel('|\rho|')
legend('|\rho_v|', '|\rho_h|',3)
grid