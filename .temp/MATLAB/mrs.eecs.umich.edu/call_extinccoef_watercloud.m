%% Source: Microwave Radar and Radiometric Remote Sensing, http://mrs.eecs.umich.edu
%% These MATLAB-based computer codes are made available to the remote
%% sensing community with no restrictions. Users may download them and
%% use them as they see fit. The codes are intended as educational tools
%% with limited ranges of applicability, so no guarantees are attached to
%% any of the codes. 
%%
clear;

%-- sample driver code (simple case): reproduces one curve in fig 8-29

t = 0; % temperature in degree C.
f = 2:1:33; % frequency in GHz
mv = 1 ; % cloud water content (g/m3)

np = length(f);
kappawc = zeros(np,1);

for ii = 1: np
    [kappawc(ii)]=ExtincCoef_WaterCloud(t,mv,f(ii));
end

figure(1)
loglog(f, kappawc)
grid




