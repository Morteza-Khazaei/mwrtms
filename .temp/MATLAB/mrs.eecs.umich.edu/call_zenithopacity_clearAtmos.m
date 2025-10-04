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
f =0:1:350   ;% frequency in GHz


nf = length(f);

tau0 = zeros(nf,1); % create empty arrays

for ii = 1: nf
    [tau0(ii)] = ZenithAtmosOpacity_ClearAtmos(T0,P0,rhov0,f(ii));
end

tau0 = tau0 .* 4.34; % transform tau to dB/km
figure(1)
semilogy(f, tau0)
grid


