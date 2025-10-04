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
f =0:1:300   ;% frequency in GHz
theta = 0; % zenith angle (degrees)

nf = length(f);

TSKY = zeros(nf,1); % create empty arrays

for ii = 1: nf
    ii
    [TSKY(ii)] = TSKY_ClearAtmos(T0,P0,rhov0,f(ii),theta);
end
%%
figure(1)
% semilogy(f, TSKY)
plot(f, TSKY)
xlabel('Frequency (GHz)')
ylabel('T_sky (K)')
grid


