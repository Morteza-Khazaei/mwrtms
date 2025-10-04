%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

T0 = 20 ;% Sea-surface temperature in 0C
P0 = 1013 ;% sea-surface total barometric pressure in mbar
rhov0= 10   ;% sea-surface water vapor density in g/m^3
f =0:1:50   ;% frequency in GHz. Max is 50 GHz
z1 =0   ;% base height of atmospheric layer in km
z2 =32   ;% top height of atmospheric layer in km
z3 = 5; % base height of the cloud layer in km
z4 = 7; % top height of the cloud layer in km
theta = 0 ;% zenith angle in degrees. Max is 80 degress
% cloud_type = 'water'; % type: 'water' or 'ice'
cloud_type = 'ice'; % type: 'water' or 'ice'
mv = 3; % cloud waer content in g/m^3. Max is 5 gm/m3


nf = length(f);

tau = zeros(nf,1); % create empty arrays
Tr = tau; % transmissivity

for ii = 1: nf
    [tau(ii) Tr(ii)] = OpticalThickn_CloudyAtmos(T0, P0, rhov0, f(ii), z1, z2, z3, z4, theta, cloud_type, mv) ;
end

tau = tau .* 4.34; % transform tau to dB/km
figure(1)
semilogy(f, tau)
xlabel('freq in GHz')
ylabel('\tau (dB/Km)')
grid

figure(2)
plot(f, Tr)
grid

