%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

f = 16; % frequency (GHz)

rho_s = 0.24; %snow density (g/cm3)

mv =  0:.1:12; % temperature ( degree C)

ri = 0.5e-3; % ice-particle radius (m)

nr = length(mv);

for n = 1:nr
    [kappa_a(n) kappa_s(n) kappa_e(n) a(n)] = MieExtinc_WetSnow(rho_s,ri,f,mv(n));
end

figure(1)

semilogy(mv, kappa_a, '-r', mv, kappa_s, '-b', mv, kappa_e, '-k')
xlabel('Liquid Water Content m_v (%)')
ylabel('\kappa (Np/m)')
legend('absorption', 'scattering', 'extinction',2)
grid
% axis([ 0 40 1.0e-4 10])


figure(2) 
semilogy(mv, a)
xlabel('Liquid Water Content m_v (%)')
ylabel('albedo')
grid