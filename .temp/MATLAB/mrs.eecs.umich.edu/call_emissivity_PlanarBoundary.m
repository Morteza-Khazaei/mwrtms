%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

eps1r = 1;

 eps2 = 54.4 -1i * 36.8;
% eps2 = 3.5 -1i * 0.2;
% eps2 = 17.9 -1i * 7.2;

theta = 0:1:90;

np = length(theta);

ev = zeros(np,1); eh = ev;

for ii = 1: np
    [eh(ii)  ev(ii)] = emissivity_PlanarBoundary(eps1r, eps2, theta(ii));
end

figure(1)

plot(theta, eh, theta, ev)
grid

%---------------------------------------

for ii = 1:np
    [rhoh(ii) rhov(ii) gammah(ii) gammav(ii) tauh(ii) tauv(ii) Th(ii) ... 
        Tv(ii)] = ReflTransm_PlanarBoundary(eps1r, eps2, theta(ii));
end

figure(2)

plot(theta, abs(rhov), theta, abs(rhoh))
xlabel('Incidence Angle (\theta_1)')
ylabel('|\rho|')
legend('|\rho_v|', '|\rho_h|',3)
grid

figure(3)

plot(theta, gammah, theta, gammav, theta, 1-gammah, theta, 1- gammav)
grid