%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 0.1:0.1:40; % frequency (GHz)

rho_s = 0.24; %snow density (g/cm3)

t = -1.0; % temperature ( degree C)

ri = 0.5e-3; % ice-particle radius (m)

nr = length(f);

for ii = 1:nr
    [kappa_a(ii) kappa_s(ii) kappa_e(ii) a(ii)] = MieExtinc_DrySnow(rho_s,ri,f(ii),t);
end

figure(1)

semilogy(f, kappa_a, '-r', f, kappa_s, '-b', f, kappa_e, '-k')
xlabel('Frequency (GHz)')
ylabel('\kappa (Np/m)')
legend('absorption', 'scattering', 'extinction',2)
grid
axis([ 0 40 1.0e-4 10])

